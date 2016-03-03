/****************************************************************************************\
 *                                                                                       *
 * Kameleon Converter                                                                    *
 *                                                                                       *
 * Copyright(c) 2016, Alexander Bock                                                     *
 * All rights reserved.                                                                  *
 *                                                                                       *
 * Redistribution and use in source and binary forms, with or without modification, are  *
 * permitted provided that the following conditions are met:                             *
 *                                                                                       *
 * 1. Redistributions of source code must retain the above copyright notice, this list   *
 * of conditions and the following disclaimer.                                           *
 *                                                                                       *
 * 2. Redistributions in binary form must reproduce the above copyright notice, this     *
 * list of conditions and the following disclaimer in the documentation and / or other   *
 * materials provided with the distribution.                                             *
 *                                                                                       *
 * 3. Neither the name of the copyright holder nor the names of its contributors may be  *
 * used to endorse or promote products derived from this software without specific prior *
 * written permission.                                                                   *
 *                                                                                       *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY   *
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES  *
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.IN NO EVENT    *
 * SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,        *
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES(INCLUDING, BUT NOT LIMITED   *
 * TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR    *
 * BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN      *
 * CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY *
 * WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH       *
 * DAMAGE.                                                                               *
 \***************************************************************************************/
 
#include <ghoul/ghoul.h>
#include <ghoul/filesystem/file.h>
#include <ghoul/filesystem/filesystem.h>
#include <ghoul/logging/logmanager.h>
#include <ghoul/logging/consolelog.h>
#include <ghoul/logging/textlog.h>
#include <ghoul/lua/ghoul_lua.h>
#include <ghoul/misc/exception.h>
#include <ccmc/Attribute.h>
#include <ccmc/Kameleon.h>

#include <algorithm>
#include <format.h>
#include <fstream>
#include <map>
#include <string>

#define _USE_MATH_DEFINES
#include <math.h>

namespace {
    const std::string KeyUseNormalization = "UseNormalization";
    bool UseNormalization = false;
}

enum class Command {
    Unknown, // This has to be the first so that the std::map selects this if the command was not found
    Native,
    Convert,
    Info
};

std::map<std::string, Command> Commands = {
    { "native", Command::Native },
    { "convert", Command::Convert },
    { "info", Command::Info }
};

enum class CoordinateSystem {
    Cartesian,
    Spherical,
    Unknown
};

namespace std {
    std::string to_string(CoordinateSystem coordinateSystem) {
        switch (coordinateSystem) {
        case CoordinateSystem::Cartesian:
            return "Cartesian";
        case CoordinateSystem::Spherical:
            return "Spherical";
        case CoordinateSystem::Unknown:
            return "Unknown";
        }
    }
} // namespace std

using namespace ghoul;

// Initializes static members
void initialize() {
    ghoul::initialize();

    logging::LogManager::initialize();
    LogMgr.addLog(std::make_shared<logging::ConsoleLog>());

    LogMgr.addLog(std::make_shared<logging::TextLog>(
        "information.txt",      // output filename
        false,                  // appending
        false,                  // time stamping
        false,                  // date stamping
        true,                   // category stamping
        false                   // logLevel stamping
    ));
}

// Deinitializes static members
void deinitialize() {
    logging::LogManager::deinitialize();

    ghoul::deinitialize();
}

void displayHelp() {
    LINFOC("Kameleon Converter", "The application can be used in one of two ways: "
        "1: By providing a CDF file as the only argument, information about this CDF "
        "file is printed to the console and stored in the 'information.txt' file."
        "2: By providing a CDF file as the first argument and a variable name contained"
        "in the file as the second argument, that variable will be extracted from the CDF"
        "file"
    );
}

void loadConfiguration(const std::string& file) {
    if (FileSys.fileExists(file)) {
        ghoul::Dictionary settings;
        ghoul::lua::loadDictionaryFromFile(file, settings);

        if (settings.hasKeyAndValue<bool>(KeyUseNormalization)) {
            UseNormalization = settings.value<bool>(KeyUseNormalization);
            LINFOC("Configuration", "Setting Normalization to " << UseNormalization);
        }
    }
}

// Returns true iff 'file' is a cdf file as determined by extension
bool isCDFFile(const std::string& file) {
    return filesystem::File(file).fileExtension() == "cdf";
}

// Initialize Kameleon with the passed 'file'
std::unique_ptr<ccmc::Kameleon> initializeKameleon(const std::string& file) {
    auto kameleon = std::make_unique<ccmc::Kameleon>();
    auto status = kameleon->open(file);
    if (status != ccmc::FileReader::OK) {
        throw RuntimeError(
            fmt::format("Error loading file {}. Status: {}", file, status),
            "Initialize Kameleon"
        );
    }
    return kameleon;
}

CoordinateSystem detectNativeCoordinateSystem(ccmc::Kameleon* kameleon) {
    bool cartesianCoordinateSystem =
        kameleon->doesVariableExist("x") &&
        kameleon->doesVariableExist("y") &&
        kameleon->doesVariableExist("z");

    bool sphericalCoordinateSystem =
        kameleon->doesVariableExist("r") &&
        kameleon->doesVariableExist("phi") &&
        kameleon->doesVariableExist("theta");

    // Only one or the other should be true
    if (!(cartesianCoordinateSystem ^ sphericalCoordinateSystem))
        throw RuntimeError(
            "Neither or both spherical or Cartesian coordinate systems are present in "
            "the file", "Kameleon Reader");
    else {
        if (cartesianCoordinateSystem)
            return CoordinateSystem::Cartesian;
        else
            return CoordinateSystem::Spherical;
    }
}

// Prints information about the passed CDF file
void printFileInformation(const std::string& file) {
    auto kameleon = initializeKameleon(file);

    // Generic information
    LINFOC("Model Name", kameleon->getModelName());
    if (kameleon->doesAttributeExist("model_type")) {
        LINFOC(
            "Model Type",
            kameleon->getGlobalAttribute("model_type").getAttributeString()
        );
    }
    LINFOC(
        "Native Coordinate System",
        std::to_string(detectNativeCoordinateSystem(kameleon.get()))
    );
    LINFOC("Time", kameleon->getCurrentTime().toString());
    LINFOC("", "");

    // Global attributes
    int nGlobalAttributes = kameleon->getNumberOfGlobalAttributes();
    LINFOC("Number of Global Attributes", nGlobalAttributes);
    for (int i = 0; i < nGlobalAttributes; ++i) {
        LINFOC(
            fmt::format("Global Attribute {}", i),
            kameleon->getGlobalAttribute(i).getAttributeName()
            );
    }
    LINFOC("", "");

    // Variable attributes
    int nVariableAttributes = kameleon->getNumberOfVariableAttributes();
    LINFOC("Number of Variable Attributes", nVariableAttributes);
    for (int i = 0; i < nVariableAttributes; ++i) {
        LINFOC(
            fmt::format("Variable Attribute {}", i),
            kameleon->getVariableAttributeName(i)
            );
    }
    LINFOC("", "");

    // Variables
    int nVariables = kameleon->getNumberOfVariables();
    LINFOC("Number of Variables", nVariables);
    LINFOC("Variables", "Name (native unit (SI unit) [min-max #nValues]");
    for (int i = 0; i < nVariables; ++i) {
        const std::string variableName = kameleon->getVariableName(i);
        const std::string siUnit = kameleon->getSIUnit(variableName);
        const std::string nativeUnit = kameleon->getNativeUnit(variableName);

        std::vector<float>* values = kameleon->getVariable(variableName);
        auto minMax = std::minmax_element(values->begin(), values->end());

        LINFOC(
            fmt::format("Variable {}", i),
            fmt::format(
                "{} ({}) ({}) [{}-{} #{}]",
                variableName,
                nativeUnit,
                siUnit,
                *minMax.first,
                *minMax.second,
                values->size()
            )
        );

        delete values;
    }
    LINFOC("", "");

    // Grids
    LINFOC(
        "Standard Grid System",
        kameleon->getGlobalAttribute("standard_grid_target").getAttributeString()
    );

    int nGrids = kameleon->getGlobalAttribute("grid_system_count").getAttributeInt();
    LINFOC("Number of Coordinate Grids", nGrids);
    for (int i = 0; i < nGrids; ++i) {
        std::string gridName = fmt::format("grid_system_{}", i + 1); // 1-based indexing

        std::string gridInformation =
            kameleon->getGlobalAttribute(gridName).getAttributeString();

        // Assuming three-dimensional data here ---abock
        const std::string xDimensionKey = gridName + "_dimension_1_size";
        int xDimension = kameleon->getGlobalAttribute(xDimensionKey).getAttributeInt();

        const std::string yDimensionKey = gridName + "_dimension_2_size";
        int yDimension = kameleon->getGlobalAttribute(yDimensionKey).getAttributeInt();

        const std::string zDimensionKey = gridName + "_dimension_3_size";
        int zDimension = kameleon->getGlobalAttribute(zDimensionKey).getAttributeInt();

        std::string dimensionInformation = fmt::format(
            "({}, {}, {})",
            xDimension,
            yDimension,
            zDimension
        );

        LINFOC(
            fmt::format("Grid {}", i),
            fmt::format("{} {}", gridInformation, dimensionInformation)
        );
    }
}

std::pair<std::string, std::string> generateFilenames(
    const std::string& file, const std::string& variable)
{
    std::string base = fmt::format("{}_{}.", filesystem::File(file).baseName(), variable);
    return std::make_pair(
        base + ".raw",
        base + ".dat"
    );
}

void convertVolume(const std::string& file, int resolution, const std::string& variable) {
    auto files = generateFilenames(file, variable);
    const std::string& outputRawFilename = files.first;
    const std::string& outputDatFilename = files.second;

    LINFOC("File Conversion",
        fmt::format("Creating files '{}' and '{}' from '{}' using variable '{}'",
            outputRawFilename,
            outputDatFilename,
            file,
            variable
            )
        );

    auto kameleon = initializeKameleon(file);

    if (!kameleon->doesVariableExist(variable)) {
        throw RuntimeError(
            fmt::format("CDF file '{}' does not contain variable '{}'", file, variable),
            "File Conversion"
            );
    }

    ccmc::Interpolator* interpolator = kameleon->createNewInterpolator();
    //ccmc::Interpolator* interpolator = kameleon->model->createNewInterpolator();

    CoordinateSystem coordinateSystem = detectNativeCoordinateSystem(kameleon.get());

    if (coordinateSystem == CoordinateSystem::Spherical) {
        // Using a Cartesian coordinate system

        // Find the minimum/maximum values of the radius
        float rMin = kameleon->getVariableAttribute("r", "actual_min").getAttributeFloat();
        float rMax = kameleon->getVariableAttribute("r", "actual_max").getAttributeFloat();
        float phiMin = kameleon->getVariableAttribute("phi", "actual_min").getAttributeFloat();
        float phiMax = kameleon->getVariableAttribute("phi", "actual_max").getAttributeFloat();
        float thetaMin = kameleon->getVariableAttribute("theta", "actual_min").getAttributeFloat();
        float thetaMax = kameleon->getVariableAttribute("theta", "actual_max").getAttributeFloat();

        LINFOC("File Conversion", "Converting values");
        std::vector<float> values(resolution * resolution * resolution);
        for (int x = 0; x < resolution; ++x) {
            for (int y = 0; y < resolution; ++y) {
                for (int z = 0; z < resolution; ++z) {
                    glm::vec3 voxelPos = glm::vec3(
                        static_cast<float>(x) / static_cast<float>(resolution - 1),
                        static_cast<float>(y) / static_cast<float>(resolution - 1),
                        static_cast<float>(z) / static_cast<float>(resolution - 1)
                    );

                    glm::vec3 pos = rMax * (2.f * voxelPos - 1.f);

                    float r = glm::length(pos);
                    float theta = acos(pos.z / r);
                    float phi = atan2(pos.y, pos.x);
                    if (phi < 0)
                        phi += M_PI * 2.f;

                    float value = 0.f;
                    if (r >= rMin && r <= rMax && 
                        theta >= thetaMin && theta <= thetaMax /*&& 
                        phi >= phiMin && phi <= phiMax*/)
                    {
                        r = r / ccmc::constants::AU_in_meters;
                        theta = (glm::degrees(theta) - 90) * -1;
                        phi = glm::degrees(phi);
                        value = interpolator->interpolate(variable, r, theta, phi);
                    }

                    int index = z * resolution * resolution + y * resolution + x;
                    values[index] = value;
                }
            }
        }

        std::vector<float>* allValues = kameleon->getVariable(variable);
        auto minMax = std::minmax_element(allValues->begin(), allValues->end());
        if (UseNormalization) {
            for (int i = 0; i < values.size(); ++i)
                values[i] = (values[i] - *minMax.first) / (*minMax.second - *minMax.first);
        }

        LINFOC("File Conversion", "Saving file '" << outputRawFilename << "'");
        std::ofstream rawOutput(outputRawFilename, std::ios_base::binary);
        rawOutput.write(
            reinterpret_cast<const char*>(values.data()),
            values.size() * sizeof(float)
            );

        LINFOC("File Conversion", "Saving file '" << outputDatFilename << "'");
        std::ofstream datOutput(outputDatFilename);
        datOutput << "RawFile: " << outputRawFilename << std::endl <<
            "Resolution: " << resolution << " " << resolution << " " << resolution << std::endl <<
            "Format: FLOAT32" << std::endl <<
            "DataRange: " << *minMax.first << " " << *minMax.second;

        delete interpolator;

    }
    else {
        throw RuntimeError("This function has not been implemented", "File Conversion");
    }

}

void extractVolume(const std::string& file, const std::string& variable) {
    auto files = generateFilenames(file, variable);
    const std::string& outputRawFilename = files.first;
    const std::string& outputDatFilename = files.second;
    
    LINFOC("File Conversion",
        fmt::format("Creating files '{}' and '{}' from '{}' using variable '{}'",
            outputRawFilename,
            outputDatFilename,
            file,
            variable
        )
    );

    auto kameleon = initializeKameleon(file);

    if (!kameleon->doesVariableExist(variable)) {
        throw RuntimeError(
            fmt::format("CDF file '{}' does not contain variable '{}'", file, variable),
            "File Conversion"
        );
    }

    ccmc::Interpolator* interpolator = kameleon->model->createNewInterpolator();
    //ccmc::Interpolator* interpolator = kameleon->createNewInterpolator();

    CoordinateSystem coordinateSystem = detectNativeCoordinateSystem(kameleon.get());

    std::string xAxis;
    std::string yAxis;
    std::string zAxis;
    switch (coordinateSystem) {
    case CoordinateSystem::Cartesian:
        xAxis = "x";
        yAxis = "y";
        zAxis = "z";
        break;
    case CoordinateSystem::Spherical:
        xAxis = "r";
        yAxis = "theta";
        zAxis = "phi";
        break;
    case CoordinateSystem::Unknown:
        throw RuntimeError("Unknown coordinate system", "File Conversion");
    }

    std::vector<float>* xAxisValues = kameleon->getVariable(xAxis);
    int xDimension = xAxisValues->size();
    std::vector<float>* yAxisValues = kameleon->getVariable(yAxis);
    int yDimension = yAxisValues->size();
    std::vector<float>* zAxisValues = kameleon->getVariable(zAxis);
    int zDimension = zAxisValues->size();
    int nValues = xDimension * yDimension * zDimension;


    LINFOC("File Conversion", "Converting values");
    std::vector<float> values(nValues);
    for (int x = 0; x < xAxisValues->size(); ++x) {
        for (int y = 0; y < yAxisValues->size(); ++y) {
            for (int z = 0; z < zAxisValues->size(); ++z) {
                // Linearize index
                int index = z * xDimension * yDimension + y * xDimension + x;
                float xPosition = xAxisValues->at(x);
                float yPosition = yAxisValues->at(y);
                float zPosition = zAxisValues->at(z);

                if (coordinateSystem == CoordinateSystem::Spherical) {
                    yPosition = (glm::degrees(yPosition)-90) * -1;
                    zPosition = glm::degrees(zPosition);
                }

                float value = interpolator->interpolate(
                    variable,
                    xPosition,
                    yPosition,
                    zPosition
                );

                values[index] = value;
            }
        }
    }

    std::vector<float>* allValues = kameleon->getVariable(variable);
    auto minMax = std::minmax_element(allValues->begin(), allValues->end());
    if (UseNormalization) {
        for (int i = 0; i < values.size(); ++i)
            values[i] = (values[i] - *minMax.first) / (*minMax.second - *minMax.first);
    }

    LINFOC("File Conversion", "Saving file '" << outputRawFilename << "'");
    std::ofstream rawOutput(outputRawFilename, std::ios_base::binary);
    rawOutput.write(
        reinterpret_cast<const char*>(values.data()),
        values.size() * sizeof(float)
    );

    LINFOC("File Conversion", "Saving file '" << outputDatFilename << "'");
    std::ofstream datOutput(outputDatFilename);
    datOutput << "RawFile: " << outputRawFilename << std::endl <<
    "Resolution: " << xDimension << " " << yDimension << " " << zDimension << std::endl <<
    "Format: FLOAT32" << std::endl <<
    "DataRange: " << *minMax.first << " " << *minMax.second;

    delete interpolator;
}

int main(int argc, char** argv) {
    const std::string ConfigurationFile = "config.cfg";

    try {
        ::initialize();

        loadConfiguration(ConfigurationFile);

        if (argc == 1)
            displayHelp();
        else {
            std::string c = argv[1];
            std::string file = argv[2];
            if (!isCDFFile(file))
                throw RuntimeError("Application must be pointed to a CDF file");

            Command command = Commands[c];

            switch (command) {
                case Command::Convert:
                {
                    std::string resolution = argv[3];

                    for (int i = 4; i < argc; ++i) {
                        std::string variable = argv[i];
                        convertVolume(file, atoi(resolution.c_str()), variable);
                    }
                    break;
                }
                case Command::Info:
                {
                    printFileInformation(file);
                    break;
                }
                case Command::Native:
                    for (int i = 3; i < argc; ++i) {
                        std::string variable = argv[i];
                        extractVolume(file, variable);
                    }
                    break;
                case Command::Unknown:
                    displayHelp();
                    break;
            }
        }

        ::deinitialize();
        std::exit(EXIT_SUCCESS);
    }
    catch (const RuntimeError& error) {
        LFATALC(error.component, error.message);
    }
}