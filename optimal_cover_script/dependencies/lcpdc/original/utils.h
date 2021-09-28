
//All this written by Simon Puglisi
//July 2006

#ifndef UTILS_H
#define UTILS_H

double getTime();

void pretty_putchar(int c);

int logbase2(unsigned int v);

void printCharArray(char *name, char *a, int n);

void printBoolArray(char *name, bool *a, int n);

void printIntArray(char *name, int *a, int n);

bool isNLE(unsigned char *x, int i, int p);

int findRunEnd(unsigned char *x, int offset, int period);

void printRun(unsigned char *x, int i, int j, int p);

bool hasPeriod1(unsigned char *x, int n);

bool match(unsigned char *a, unsigned char *b, int n);

int longestmatch(unsigned char *x, int i, int j, int n);

#endif

