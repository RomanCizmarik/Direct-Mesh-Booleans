

#ifdef _MSC_VER // Workaround for known bugs and issues on MSVC
#define _HAS_STD_BYTE 0  // https://developercommunity.visualstudio.com/t/error-c2872-byte-ambiguous-symbol/93889
#define NOMINMAX // https://stackoverflow.com/questions/1825904/error-c2589-on-stdnumeric-limitsdoublemin
#endif

#include <iostream>

#include "SimpleMesh.h"

#include "Booleans.h"

void FWNTreeTest(DMB::SimpleMesh& m)
{
    //generate random points
    std::vector<DMB::SimpleMesh::Point> randomPoints;
    {
        Eigen::MatrixXd V;
        Eigen::MatrixXi F;

        DMB::OpenMesh2Matrix(m, V, F);

        int nPoints = 100;

        // Generate a list of random query points in the bounding box
        Eigen::MatrixXd Q = Eigen::MatrixXd::Random(nPoints, 3);
        const Eigen::RowVector3d Vmin = V.colwise().minCoeff();
        const Eigen::RowVector3d Vmax = V.colwise().maxCoeff();
        const Eigen::RowVector3d Vdiag = Vmax - Vmin;
        for (int q = 0; q < Q.rows(); q++)
        {
            Q.row(q) = (Q.row(q).array() * 0.5 + 0.5) * Vdiag.array() + Vmin.array();
        }

        for (int q = 0; q < Q.rows(); q++)
        {
            randomPoints.push_back({ Q.row(q).coeff(0),Q.row(q).coeff(1),Q.row(q).coeff(2) });
        }
    }

    DMB::FastWindingNumber< DMB::SimpleMesh > tree;

    m.update_face_normals();
    tree.build(m);

    for (auto i = 0; i < randomPoints.size(); ++i)
    {
        auto fastResult = tree.windingNumber(randomPoints[i]);
        auto exactResult = DMB::windingNumber(m, randomPoints[i]);

        auto insideFast = fastResult > 0.5;
        auto insideExact = exactResult > 0.5;

        //EXPECT_NEAR(fastResult, exactResult, 0.045);


        //if (std::abs(fastResult - exactResult) > 1e-2)
        if(insideFast != insideExact)
        {
            std::cout << "FWN approximation incorrect!" << std::endl;
            //std::cout << "fast wn: " << fastResult << std::endl;
            //std::cout << "wn: " << exactResult << std::endl;
            std::cout << "difference: " << std::abs(fastResult - exactResult) << std::endl;
        }
    }

}

template<typename FunctionType>
void time(std::string taskName, const FunctionType& func)
{
    auto t_before = clock();
    func();
    const double t_after = clock();

    std::cout << "Task " + taskName + " took: " << ((t_after - t_before) / (CLOCKS_PER_SEC / 1000)) << " miliseconds" << std::endl;
};

void FWNSpeedTest(DMB::SimpleMesh& m)
{
    //generate random points
    std::vector<DMB::SimpleMesh::Point> randomPoints;
    {
        Eigen::MatrixXd V;
        Eigen::MatrixXi F;

        DMB::OpenMesh2Matrix(m, V, F);

        int nPoints = 100;

        // Generate a list of random query points in the bounding box
        Eigen::MatrixXd Q = Eigen::MatrixXd::Random(nPoints, 3);
        const Eigen::RowVector3d Vmin = V.colwise().minCoeff();
        const Eigen::RowVector3d Vmax = V.colwise().maxCoeff();
        const Eigen::RowVector3d Vdiag = Vmax - Vmin;
        for (int q = 0; q < Q.rows(); q++)
        {
            Q.row(q) = (Q.row(q).array() * 0.5 + 0.5) * Vdiag.array() + Vmin.array();
        }

        for (int q = 0; q < Q.rows(); q++)
        {
            randomPoints.push_back({ Q.row(q).coeff(0),Q.row(q).coeff(1),Q.row(q).coeff(2) });
        }
    }

    DMB::FastWindingNumber< DMB::SimpleMesh > tree;


    time("FWN init", [&]
        {
            m.update_face_normals();
            tree.build(m);
        });

    time("FWN compute", [&]
        {
            for (auto i = 0; i < randomPoints.size(); ++i)
            {
                auto fastResult = tree.windingNumber(randomPoints[i]);
            }
        });

    time("Exact WN compute", [&]
        {
            for (auto i = 0; i < randomPoints.size(); ++i)
            {
                auto exactResult = DMB::windingNumber(m, randomPoints[i]);
            }
        });
}


//int main(int argc, char** argv)
//{
//    //for (auto fn : { "C:/skola/PhD/Samples/booleans/sphere_9x6.obj" , "C:/skola/PhD/VUT/booleans_paper/evaluation/Validity/input_models/jaw.obj" })
//    //{
//    //    std::cout << fn << std::endl;
//
//    //    DMB::SimpleMesh m;
//    //    OpenMesh::IO::read_mesh(m, fn);
//    //    FWNTreeTest(m);
//    //    FWNSpeedTest(m);
//    //}
//
//    //DMB::SimpleMesh m;
//    //DMB::meshUnion(m, "C:/skola/PhD/VUT/booleans_paper/evaluation/Validity/input_models/101954.stl", "C:/skola/PhD/VUT/booleans_paper/evaluation/Validity/input_models/894835.stl");
//    //OpenMesh::IO::write_mesh(m, "out.obj");
//
//
//
//    return 0;
//}


bool printHelp(const std::map<std::string, std::string>& arguments)
{
    if (arguments.find("h") != arguments.end() || arguments.find("help") != arguments.end())
    {
        std::cout <<
            "This executable is ment for running an algorithm as a separate process. \n"
            "To run an algorithm, share your data via shared memory and specify the algorithm to run using -algorithm parameter. \n"
            "\n"
            "Parameters: \n"
            "-algorithm:NAME_OF_ALGORITHM_TO_RUN    Runs algorithm specified by the string NAME_OF_ALGORITHM_TO_RUN. \n"
            "-h                                     Shows this help. \n"
            "-help                                  Shows this help. \n"
            "Example usage: \n"
            "./BSPWorker.exe -algorithm:BooleanUnion"
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

