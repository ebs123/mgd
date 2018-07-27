#if defined(WIN32) | defined(_WIN32) | defined(WIN64) | defined(_WIN64) 
#include "includes\NTracer.h"
#else
#include "includes/NTracer.h"
#endif

FILE *NTracer::file;

NTracer::NTracer(void)
{
}

NTracer::~NTracer(void)
{
}