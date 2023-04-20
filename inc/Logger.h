#ifndef LOGGER_H_
#define LOGGER_H_

#include "common.h"

void LogAll(const String mylog, SharedPtr<CommGrid> commgrid);
String ProcessorName(SharedPtr<CommGrid> commgrid);

#endif
