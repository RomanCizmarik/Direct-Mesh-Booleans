//TODO: License

//GTest
#include <gtest/gtest.h>

//Eigen
#include "Eigen/Dense"

//Local
#include "SimpleMesh.h"
#include "Utils.h"
#include "FastWindingNumber.h"

TEST(FWNTest, PrecisionTest)
{
    std::vector<std::string> filenamesList;

    for (auto fn : filenamesList)
    {
        DMB::SimpleMesh m;
        OpenMesh::IO::read_mesh(m, fn);
        // 
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
            if (insideFast != insideExact)
            {
                std::cout << "FWN approximation incorrect!" << std::endl;
                //std::cout << "fast wn: " << fastResult << std::endl;
                //std::cout << "wn: " << exactResult << std::endl;
                std::cout << "difference: " << std::abs(fastResult - exactResult) << std::endl;
            }
        }
    }

    
}
