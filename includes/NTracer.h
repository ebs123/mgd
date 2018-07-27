#pragma once
#include <iostream>

class NTracer
{
private:
	static FILE *file;

public:
	NTracer(void);

	virtual ~NTracer(void);

	static void traceToTerminal(const char *methodName, const int line, const char *text);
	template <class T> static void traceToFile(T *data, const char *dataType, const int &dataArraySize);
	static int openFile(const char *fileName);
	static void closeFile();
};

template <class T> void NTracer::traceToFile(T *data, const char *dataType, const int &dataArraySize)
{
	int counter = 0;

	if(!strcmp(dataType, "int"))
	{
		while(counter < dataArraySize)
		{
			fprintf(file, "%d ", *data);
			counter++;
			data++;
		}
		fprintf(file, "\n");
	}
	else if(!strcmp(dataType, "double"))
	{
		while(counter < dataArraySize)
		{
			fprintf(file, "%lf ", *data);
			counter++;
			data++;
		}
		fprintf(file, "\n");
	}

	for(counter; counter > 0; counter--)
		data--;

	return;
}

inline void NTracer::traceToTerminal(const char *methodName, const int line, const char *text)
{
	std::cout << "in method " << methodName << ", from line " << line << ", we get the error: " << text << "\n";
}

inline int NTracer::openFile(const char *fileName)
{
	file = fopen(fileName, "w");

	if(file)
		return 1;
	else
		return 0;
}

inline void NTracer::closeFile()
{
	fclose(file);
}
