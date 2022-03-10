<?php

namespace app\modules\bundle\controllers;

use common\models\billing\user\limits\Limits;
use common\models\Bundle;
use common\models\bundle\Exception\StatusException;
use common\models\EventLog;
use common\models\PartnerException;
use common\models\User;
use common\utils\ArrayAction;
use Exception;
use frontend\helpers\BundleStatus;
use frontend\models\Controller;
use frontend\models\Layout;
use frontend\models\ModalDialog\Buttons\Close;
use frontend\models\ModalDialog\Response;
use Yii;
use yii\web\ForbiddenHttpException;

class ActionController extends Controller
{
    public function actionModalAjaxPlay($bundleIds)
    {
        $errors  = [];
        $message = '';
        $started = false;

        $response        = new Response();
        $response->title = Yii::t('front', 'Run bundle');
        $closeButton     = new Close();

        if ($bundleIds !== '')
        {
            $bundleIds = explode(',', $bundleIds);
        }
        else
        {
            $response->body   = $this->renderAjax('action', [
                'started' => $started,
                'errors'  => [Yii::t('front', 'Error') => Yii::t('front', 'No bundle selected')],
                'message' => false,
            ]);
            $response->footer = [$closeButton];

            return $response;
        }

        /** @var Bundle $bundle */
        $bundles = Bundle::findAll(['id' => $bundleIds, 'userId' => Yii::$app->user->getRealId(), 'dateDeleted' => null]);

        if (is_null($bundles))
        {
            throw new ForbiddenHttpException();
        }

        /** @var User $user */
        $user           = Yii::$app->user->identity;
        $userLimits     = new Limits($user->id);
        $successBundles = [];
        $failedBundles  = [];

        foreach ($bundles as $bundle)
        {
            if (! $userLimits->canStartBundles($bundle))
            {
                if (Layout::isWhiteLabel())
                {
                    $errors['Plan'] = Yii::t('front', 'Sorry, bundle start is not allowed due to service unavailability.');
                }
                else
                {
                    $errors['Plan'] = Yii::t('front', 'Sorry, your current plan does not allow bundle start.');
                }
            }
            else
            {
                try
                {
                    $trigger       = $bundle->getTrigger();
                    $isHistoricRun = $trigger && ! in_array(Bundle::LAUNCH_MODE_NORMAL, $trigger->getBundleLaunchModes()) && $user->isAdmin();

                    if ($isHistoricRun)
                    {
                        $modes = $trigger->historicRunModes;

                        $bundle->launchMode         = Bundle::LAUNCH_MODE_HISTORIC;
                        $bundle->historicRunMode    = ArrayAction::arrayKeyFirst($modes);
                        $bundle->historicRunSubMode = array_shift($modes[$bundle->historicRunMode]);
                    }

                    $bundle->play();

                    /** @var EventLog $historyLog */
                    $historyLog = Yii::$container->get('EventLog');
                    $historyLog->setBundleId($bundle->id);
                    $historyLog->playBundle();

                    $started                     = true;
                    $successBundles[$bundle->id] = [
                        'statusStr'    => BundleStatus::widget($bundle->status, $bundle->blocked),
                        'actionTitle'  => Yii::t('front', 'Pause bundle'),
                        'statusAction' => 'pause',
                    ];
                }
                catch (StatusException | PartnerException $e)
                {
                    Yii::error($e->__toString());

                    $errors[$bundle->getReadableTitle()] = $e->getMessage();
                }
                catch (Exception $e)
                {
                    Yii::error($e->__toString());

                    $errors['Unknown'] = Yii::t('common', 'Unknown error');
                }
            }
        }
        $successBundlesCount = count($successBundles);

        if ($successBundlesCount)
        {
            $message = Yii::t('front', "{0,plural,=1{1 bundle} other{# bundles}} started", [$successBundlesCount]);

            $closeButton->setOnClickEvent([
                'event'          => 'refresh',
                'successBundles' => $successBundles,
                'failedBundles'  => $failedBundles,
            ]);
        }

        $response->body   = $this->renderAjax('action', [
            'errors'  => $errors,
            'message' => $message,
            'started' => $started,
        ]);
        $response->footer = [$closeButton];

        return $response;
    }

    public function actionModalAjaxPause($bundleIds)
    {
        $errors  = [];
        $message = '';

        $response        = new Response();
        $response->title = Yii::t('front', 'Pause bundle');

        $closeButton = new Close();

        if ($bundleIds !== '')
        {
            $bundleIds = explode(',', $bundleIds);
        }
        else
        {
            $response->body   = $this->renderAjax('action', [
                'errors'  => [Yii::t('front', 'Error') => Yii::t('front', 'No bundle selected')],
                'message' => false,
                'started' => false,
            ]);
            $response->footer = [$closeButton];

            return $response;
        }

        /** @var Bundle $bundle */
        $bundles = Bundle::findAll(['id' => $bundleIds, 'userId' => Yii::$app->user->getRealId(), 'dateDeleted' => null]);

        if (is_null($bundles))
        {
            throw new ForbiddenHttpException();
        }

        $successBundles = [];
        $failedBundles  = [];
        foreach ($bundles as $bundle)
        {
            try
            {
                $bundle->pause();

                /** @var EventLog $historyLog */
                $historyLog = Yii::$container->get('EventLog');
                $historyLog->setBundleId($bundle->id);
                $historyLog->pauseBundle();

                $successBundles[$bundle->id] = [
                    'statusStr'    => BundleStatus::widget($bundle->status, $bundle->blocked),
                    'actionTitle'  => Yii::t('front', 'Run bundle'),
                    'statusAction' => 'play',
                ];
            }
            catch (StatusException | PartnerException $e)
            {
                Yii::error($e->__toString());

                $errors[$bundle->getReadableTitle()] = $e->getMessage();
            }
            catch (Exception $e)
            {
                Yii::error($e->__toString());

                $errors['Unknown'] = Yii::t('common', 'Unknown error');
            }
        }
        $successBundlesCount = count($successBundles);

        if ($successBundlesCount)
        {
            $message = Yii::t('front', "{0,plural,=1{1 bundle} other{# bundles}} suspended", [$successBundlesCount]);

            $closeButton->setOnClickEvent([
                'event'          => 'refresh',
                'successBundles' => $successBundles,
                'failedBundles'  => $failedBundles,
            ]);
        }

        $response->footer = [$closeButton];
        $response->body   = $this->renderAjax('action', [
            'errors'  => $errors,
            'message' => $message,
            'started' => false,
        ]);

        return $response;
    }

    public function actionModalAjaxRestore($bundleIds)
    {
        $errors  = [];
        $message = '';
        $started = false;

        $response        = new Response();
        $closeButton     = new Close();
        $response->title = Yii::t('front', 'Restore deleted');
        if ($bundleIds !== '')
        {
            $bundleIds = explode(',', $bundleIds);
            $bundles   = Bundle::find()
                ->where(['id' => $bundleIds, 'userId' => Yii::$app->user->getRealId()])
                ->andWhere(['not', ['dateDeleted' => null]])
                ->all();
            foreach ($bundles as $bundle)
            {
                $bundle->dateDeleted = null;
                $bundle->save();
            }
            $message = Yii::t('front', '{0,plural,=1{1 bundle} other{# bundles}} restored', count($bundles));
            $closeButton->setOnClickEvent([
                'event'     => 'refresh',
                'bundleIds' => $bundleIds,
            ]);
        }
        $response->body   = $this->renderAjax('action', [
            'errors'  => $errors,
            'message' => $message,
            'started' => $started,

        ]);
        $response->footer = [$closeButton];

        return $response;
    }
}
