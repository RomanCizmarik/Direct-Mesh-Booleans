

#ifdef _MSC_VER // Workaround for known bugs and issues on MSVC
#define _HAS_STD_BYTE 0  // https://developercommunity.visualstudio.com/t/error-c2872-byte-ambiguous-symbol/93889
#define NOMINMAX // https://stackoverflow.com/questions/1825904/error-c2589-on-stdnumeric-limitsdoublemin
#endif

#include <iostream>

#include "SimpleMesh.h"

#include "Booleans.h"


bool printHelp(const std::map<std::string, std::string>& arguments)
{
    if (arguments.find("h") != arguments.end() || arguments.find("help") != arguments.end())
    {
        std::cout <<
            "This executable loads two input meshes and computes requested Boolean operation (Union, Difference, or Intersection). " 
            "If output file path is provided, the result will be saved to the file specified. \n"
            "Supported formats for input and output files are: .OBJ, .OFF, .STL.\n"
            "Arguments of this executable:\n"
            "[U|I|D] : Specifies the Boolean operation to compute - Union, Intersection, or Difference\n"
            "input1 - Path to the first input mesh file.\n"
            "input2 - Path to the second input mesh file.\n"
            "output - Path to the output mesh file.\n"
            "Usage:\n"
            "DirectMeshBooleansExample.exe [U|D|I] input1 input2 outputPath\n"
            << std::endl;
        return true;
    }

    return false;
}

void parseArguments(int argc, char** argv, std::map<std::string, std::string>& arguments, std::vector<std::string>& unnamedArguments)
{
    const std::string delimiter = ":";
    const std::string prefix = "-";

    for (int a = 1; a < argc; ++a)
    {
        std::string s(argv[a]);

        //unnamed arguments - not starting with prefix
        if (s.substr(0, 1) != prefix)
        {
            unnamedArguments.push_back(s);
            continue;
        }

        auto it = s.find(delimiter);
        std::string key = s.substr(1, it - 1); //get rid of prefix
        std::string value = s.substr(it + 1, s.size()); //get rid of delimiter

        arguments[key] = value;
    }
}

bool runAlgorithm(const std::map<std::string, std::string>& arguments, const std::vector<std::string>& unnamedArguments)
{
    DMB::SimpleMesh result;
    bool success = false;

    if (unnamedArguments[0] == "U")
    {
        success = DMB::meshUnion(result, unnamedArguments[1], unnamedArguments[2]);
    }
    else if (unnamedArguments[0] == "D")
    {
        success = DMB::meshSubtraction(result, unnamedArguments[1], unnamedArguments[2]);
    }
    else if (unnamedArguments[0] == "I")
    {
        success = DMB::meshIntersection(result, unnamedArguments[1], unnamedArguments[2]);
    }
    else
    {
        return false;
    }

    if (success && unnamedArguments.size() > 3)
    {
        try
        {
            OpenMesh::IO::write_mesh(result, unnamedArguments[3]);
        }
        catch (...)
        {
            std::cout << "writing failed" << std::endl;

        }
    }

    return success;
}

int main(int argc, char** argv)
{
    std::map<std::string, std::string> namedArguments;
    std::vector<std::string> unnamedArguments;

    parseArguments(argc, argv, namedArguments, unnamedArguments);

    if (printHelp(namedArguments))
    {
        return 0;
    }

    try
    {
        if (!runAlgorithm(namedArguments, unnamedArguments))
        {
            return 1;
        }

    }
    catch (std::exception& e)
    {
        std::cerr << "Fatal exception: " << e.what() << std::endl;
        return 1;
    }
    catch (...)
    {
        std::cerr << "Unknown exception" << std::endl;
        return 1;
    }

    return 0;
}

