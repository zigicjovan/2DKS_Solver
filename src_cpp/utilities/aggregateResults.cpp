#include <algorithm>
#include <array>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

using namespace std;
namespace fs = std::filesystem;

using Solution = array<double, 7>;
using Solutions = vector<Solution>;

enum class FileType {
    Branch,
    InitialEnergyPowerLaw,
    DomainSizePowerLaw,
    EnergyTimeWindowPowerLaw,
    DomainTimeWindowPowerLaw,
    Unknown
};

bool beginsWith(const string& text, const string& prefix)
{
    return text.compare(0, prefix.size(), prefix) == 0;
}

FileType getFileType(const string& filename)
{
    if (beginsWith(filename, "branch_IC_"))
        return FileType::Branch;

    if (beginsWith(filename, "powerlawK_IC_"))
        return FileType::InitialEnergyPowerLaw;

    if (beginsWith(filename, "powerlawL_IC_"))
        return FileType::DomainSizePowerLaw;

    if (beginsWith(filename, "powerlawTK_IC_"))
        return FileType::EnergyTimeWindowPowerLaw;

    if (beginsWith(filename, "powerlawTL_IC_"))
        return FileType::DomainTimeWindowPowerLaw;

    return FileType::Unknown;
}

bool solutionsMatch(
    const Solution& first,
    const Solution& second,
    FileType fileType
) {
    switch (fileType) {
        case FileType::Branch:
            return first[3] == second[3];

        case FileType::InitialEnergyPowerLaw:
        case FileType::EnergyTimeWindowPowerLaw:
            return first[0] == second[0];

        case FileType::DomainSizePowerLaw:
        case FileType::DomainTimeWindowPowerLaw:
            return first[1] == second[1] &&
                   first[2] == second[2];

        case FileType::Unknown:
            return false;
    }

    return false;
}

bool solutionComesBefore(
    const Solution& first,
    const Solution& second,
    FileType fileType
) {
    switch (fileType) {
        case FileType::Branch:
            return first[3] < second[3];

        case FileType::InitialEnergyPowerLaw:
        case FileType::EnergyTimeWindowPowerLaw:
            return first[0] < second[0];

        case FileType::DomainSizePowerLaw:
        case FileType::DomainTimeWindowPowerLaw:
            if (first[1] != second[1])
                return first[1] < second[1];

            return first[2] < second[2];

        case FileType::Unknown:
            return false;
    }

    return false;
}

void insertOrReplace(
    Solutions& solutions,
    const Solution& candidate,
    FileType fileType
) {
    auto existing = find_if(
        solutions.begin(),
        solutions.end(),
        [&](const Solution& solution) {
            return solutionsMatch(
                solution,
                candidate,
                fileType
            );
        }
    );

    if (existing == solutions.end()) {
        solutions.push_back(candidate);
    }
    else if (candidate[6] > (*existing)[6]) {
        *existing = candidate;
    }
}

void readSolutionFile(
    const fs::path& inputPath,
    Solutions& solutions,
    FileType fileType
) {
    ifstream inputFile(inputPath);

    if (!inputFile) {
        cerr << "Warning: could not open "
             << inputPath << '\n';
        return;
    }

    Solution candidate;

    while (
        inputFile
        >> candidate[0]
        >> candidate[1]
        >> candidate[2]
        >> candidate[3]
        >> candidate[4]
        >> candidate[5]
        >> candidate[6]
    ) {
        insertOrReplace(
            solutions,
            candidate,
            fileType
        );
    }

    if (!inputFile.eof()) {
        cerr << "Warning: incomplete or invalid row in "
             << inputPath << '\n';
    }
}

void writeSolutionFile(
    const fs::path& outputPath,
    Solutions& solutions,
    FileType fileType
) {
    sort(
        solutions.begin(),
        solutions.end(),
        [&](const Solution& first, const Solution& second) {
            return solutionComesBefore(
                first,
                second,
                fileType
            );
        }
    );

    ofstream outputFile(outputPath);

    if (!outputFile) {
        throw runtime_error(
            "Could not create output file: " +
            outputPath.string()
        );
    }

    outputFile << scientific << setprecision(16);

    for (const Solution& solution : solutions) {
        for (size_t column = 0;
             column < solution.size();
             ++column) {

            if (column > 0)
                outputFile << ' ';

            outputFile << solution[column];
        }

        outputFile << '\n';
    }

    if (!outputFile) {
        throw runtime_error(
            "Error while writing output file: " +
            outputPath.string()
        );
    }
}

int main(int argc, char* argv[])
{
    if (argc != 3) {
        cerr << "Usage:\n"
             << "  " << argv[0]
             << " INPUT_DIRECTORY OUTPUT_DIRECTORY\n\n"
             << "Example:\n"
             << "  " << argv[0]
             << " Data Aggregated\n";

        return 1;
    }

    const fs::path inputDirectory = argv[1];
    const fs::path outputDirectory = argv[2];

    const fs::path absoluteOutputDirectory = fs::absolute(outputDirectory).lexically_normal();

    if (!fs::exists(inputDirectory)) {
        cerr << "Input directory does not exist: "
             << inputDirectory << '\n';
        return 1;
    }

    if (!fs::is_directory(inputDirectory)) {
        cerr << "Input path is not a directory: "
             << inputDirectory << '\n';
        return 1;
    }

    map<string, Solutions> aggregatedSolutions;
    map<string, FileType> fileTypes;

    try {
        for (const fs::directory_entry& entry :
            fs::recursive_directory_iterator(inputDirectory)) {

            if (!entry.is_regular_file())
                continue;

            const fs::path absoluteInputPath =
                fs::absolute(entry.path()).lexically_normal();

            const string filename =
                entry.path().filename().string();

            const FileType fileType =
                getFileType(filename);

            if (fileType == FileType::Unknown)
                continue;

            fileTypes[filename] = fileType;

            readSolutionFile(
                entry.path(),
                aggregatedSolutions[filename],
                fileType
            );
        }

        fs::create_directories(outputDirectory);

        for (auto& filenameAndSolutions :
             aggregatedSolutions) {

            const string& filename =
                filenameAndSolutions.first;

            Solutions& solutions =
                filenameAndSolutions.second;

            const FileType fileType =
                fileTypes.at(filename);

            const fs::path outputPath =
                outputDirectory / filename;

            writeSolutionFile(
                outputPath,
                solutions,
                fileType
            );
        }
    }
    catch (const fs::filesystem_error& error) {
        cerr << "Filesystem error: "
             << error.what() << '\n';
        return 1;
    }
    catch (const exception& error) {
        cerr << "Aggregation failed: "
             << error.what() << '\n';
        return 1;
    }

    cout << "Aggregation complete: "
         << aggregatedSolutions.size()
         << " files created in "
         << outputDirectory << '\n';

    return 0;
}