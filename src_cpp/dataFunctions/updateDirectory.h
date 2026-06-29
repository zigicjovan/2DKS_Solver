#ifndef UPDATEDIRECTORYDATA_H
#define UPDATEDIRECTORYDATA_H

#include "Pathnames.h"
#include "SolutionData.h"

#include <filesystem>

std::filesystem::path appendTimeStep(const std::filesystem::path& path, double dCurrentT);
void loadData(const Pathnames& paths, SolutionData& storedData, SolutionDataType storedDataType, double dCurrentT = 0.0);
void saveData(const Pathnames& paths, SolutionData& storedData, SolutionDataType storedDataType, double dCurrentT = 0.0);
void deleteData(const Pathnames& paths);

#endif