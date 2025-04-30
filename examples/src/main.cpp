

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


int main(int argc, char** argv)
{



	////OpenMesh::IO::read_mesh(m, "C:/skola/PhD/Samples/booleans/result19.obj");
//    m.update_normals();

    

    for (auto fn : { "C:/skola/PhD/Samples/booleans/sphere_9x6.obj" , "C:/skola/PhD/VUT/booleans_paper/evaluation/Validity/input_models/jaw.obj" })
    {
        std::cout << fn << std::endl;

        DMB::SimpleMesh m;
        OpenMesh::IO::read_mesh(m, fn);
        FWNTreeTest(m);
        FWNSpeedTest(m);
    }

    //DMB::SimpleMesh m;
    //DMB::meshUnion(m, "C:/skola/PhD/VUT/booleans_paper/evaluation/Validity/input_models/101954.stl", "C:/skola/PhD/VUT/booleans_paper/evaluation/Validity/input_models/894835.stl");
    

}