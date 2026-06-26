#ifndef UPDATEDIRECTORYDATA_H
#define UPDATEDIRECTORYDATA_H

#include "Pathnames.h"
#include "SolutionData.h"

std::filesystem::path appendTimeStep(const std::filesystem::path &path, int currentT);
void loadData(const Pathnames &paths, SolutionData &storedData, SolutionDataType storedDataType, int currentT = 0);
void saveData(const Pathnames &paths, SolutionData &storedData, SolutionDataType storedDataType, int currentT = 0);
void deleteData(const Pathnames &paths);

#endif