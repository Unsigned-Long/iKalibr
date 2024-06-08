#ifndef PCL_GICP_OMP_H_
#define PCL_GICP_OMP_H_

#include "atomic"
#include "pcl/registration/bfgs.h"
#include "pcl/registration/boost.h"
#include "pcl/registration/exceptions.h"
#include "pcl/registration/icp.h"
#include "pclomp/gicp_omp.hpp"

namespace pclomp {
    /** \brief GeneralizedIterativeClosestPoint is an ICP variant that implements the
      * generalized iterative closest point algorithm as described by Alex Segal et al. in
      * http://www.robots.ox.ac.uk/~avsegal/resources/papers/Generalized_ICP.pdf
      * The approach is based on using anistropic cost functions to optimize the alignment
      * after closest point assignments have been made.
      * The original code uses GSL and ANN while in ours we use an eigen mapped BFGS and
      * FLANN.
      * \author Nizar Sallem
      * \ingroup registration
      */
    template<typename PointSource, typename PointTarget>
    class GeneralizedIterativeClosestPoint : public pcl::IterativeClosestPoint<PointSource, PointTarget> {
    public:
        using pcl::IterativeClosestPoint<PointSource, PointTarget>::reg_name_;
        using pcl::IterativeClosestPoint<PointSource, PointTarget>::getClassName;
        using pcl::IterativeClosestPoint<PointSource, PointTarget>::indices_;
        using pcl::IterativeClosestPoint<PointSource, PointTarget>::target_;
        using pcl::IterativeClosestPoint<PointSource, PointTarget>::input_;
        using pcl::IterativeClosestPoint<PointSource, PointTarget>::tree_;
        using pcl::IterativeClosestPoint<PointSource, PointTarget>::tree_reciprocal_;
        using pcl::IterativeClosestPoint<PointSource, PointTarget>::nr_iterations_;
        using pcl::IterativeClosestPoint<PointSource, PointTarget>::max_iterations_;
        using pcl::IterativeClosestPoint<PointSource, PointTarget>::previous_transformation_;
        using pcl::IterativeClosestPoint<PointSource, PointTarget>::final_transformation_;
        using pcl::IterativeClosestPoint<PointSource, PointTarget>::transformation_;
        using pcl::IterativeClosestPoint<PointSource, PointTarget>::transformation_epsilon_;
        using pcl::IterativeClosestPoint<PointSource, PointTarget>::converged_;
        using pcl::IterativeClosestPoint<PointSource, PointTarget>::corr_dist_threshold_;
        using pcl::IterativeClosestPoint<PointSource, PointTarget>::inlier_threshold_;
        using pcl::IterativeClosestPoint<PointSource, PointTarget>::min_number_correspondences_;
        using pcl::IterativeClosestPoint<PointSource, PointTarget>::update_visualizer_;

        using PointCloudSource = pcl::PointCloud<PointSource>;
        using PointCloudSourcePtr = typename PointCloudSource::Ptr;
        using PointCloudSourceConstPtr = typename PointCloudSource::ConstPtr;

        using PointCloudTarget = pcl::PointCloud<PointTarget>;
        using PointCloudTargetPtr = typename PointCloudTarget::Ptr;
        using PointCloudTargetConstPtr = typename PointCloudTarget::ConstPtr;

        using PointIndicesPtr = pcl::PointIndices::Ptr;
        using PointIndicesConstPtr = pcl::PointIndices::ConstPtr;

        using MatricesVector = std::vector<Eigen::Matrix3d, Eigen::aligned_allocator<Eigen::Matrix3d> >;
        using MatricesVectorPtr = boost::shared_ptr<MatricesVector>;
        using MatricesVectorConstPtr = boost::shared_ptr<const MatricesVector>;

        using InputKdTree = typename pcl::Registration<PointSource, PointTarget>::KdTree;
        using InputKdTreePtr = typename pcl::Registration<PointSource, PointTarget>::KdTreePtr;

        using Ptr = boost::shared_ptr<GeneralizedIterativeClosestPoint<PointSource, PointTarget> >;
        using ConstPtr = boost::shared_ptr<const GeneralizedIterativeClosestPoint<PointSource, PointTarget> >;


        using Vector6d = Eigen::Matrix<double, 6, 1>;

        /** \brief Empty constructor. */
        GeneralizedIterativeClosestPoint()
                : k_correspondences_(20), gicp_epsilon_(0.001), rotation_epsilon_(2e-3), mahalanobis_(0),
                  max_inner_iterations_(20) {
            min_number_correspondences_ = 4;
            reg_name_ = "GeneralizedIterativeClosestPoint";
            max_iterations_ = 200;
            transformation_epsilon_ = 5e-4;
            corr_dist_threshold_ = 5.;
            rigid_transformation_estimation_ = [this](const PointCloudSource &cloud_src,
                                                      const std::vector<int> &indices_src,
                                                      const PointCloudTarget &cloud_tgt,
                                                      const std::vector<int> &indices_tgt,
                                                      Eigen::Matrix4f &transformation_matrix) {
                estimateRigidTransformationBFGS(cloud_src, indices_src, cloud_tgt, indices_tgt, transformation_matrix);
            };
        }

        /** \brief Provide a pointer to the input dataset
          * \param cloud the const boost shared pointer to a PointCloud message
          */
        inline void
        setInputSource(const PointCloudSourceConstPtr &cloud) override {

            if (cloud->points.empty()) {
                PCL_ERROR ("[pcl::%s::setInputSource] Invalid or empty point cloud dataset given!\n",
                           getClassName().c_str());
                return;
            }
            PointCloudSource input = *cloud;
            // Set all the point.data[3] values to 1 to aid the rigid transformation
            for (size_t i = 0; i < input.size(); ++i)
                input[i].data[3] = 1.0;

            pcl::IterativeClosestPoint<PointSource, PointTarget>::setInputSource(cloud);
            input_covariances_.reset();
        }

        /** \brief Provide a pointer to the covariances of the input source (if computed externally!).
          * If not set, GeneralizedIterativeClosestPoint will compute the covariances itself.
          * Make sure to set the covariances AFTER setting the input source point cloud (setting the input source point cloud will reset the covariances).
          * \param[in] target the input point cloud target
          */
        inline void
        setSourceCovariances(const MatricesVectorPtr &covariances) {
            input_covariances_ = covariances;
        }

        /** \brief Provide a pointer to the input target (e.g., the point cloud that we want to align the input source to)
          * \param[in] target the input point cloud target
          */
        inline void
        setInputTarget(const PointCloudTargetConstPtr &target) override {
            pcl::IterativeClosestPoint<PointSource, PointTarget>::setInputTarget(target);
            target_covariances_.reset();
        }

        /** \brief Provide a pointer to the covariances of the input target (if computed externally!).
          * If not set, GeneralizedIterativeClosestPoint will compute the covariances itself.
          * Make sure to set the covariances AFTER setting the input source point cloud (setting the input source point cloud will reset the covariances).
          * \param[in] target the input point cloud target
          */
        inline void
        setTargetCovariances(const MatricesVectorPtr &covariances) {
            target_covariances_ = covariances;
        }

        /** \brief Estimate a rigid rotation transformation between a source and a target point cloud using an iterative
          * non-linear Levenberg-Marquardt approach.
          * \param[in] cloud_src the source point cloud dataset
          * \param[in] indices_src the vector of indices describing the points of interest in \a cloud_src
          * \param[in] cloud_tgt the target point cloud dataset
          * \param[in] indices_tgt the vector of indices describing the correspondences of the interest points from \a indices_src
          * \param[out] transformation_matrix the resultant transformation matrix
          */
        void
        estimateRigidTransformationBFGS(const PointCloudSource &cloud_src,
                                        const std::vector<int> &indices_src,
                                        const PointCloudTarget &cloud_tgt,
                                        const std::vector<int> &indices_tgt,
                                        Eigen::Matrix4f &transformation_matrix);

        /** \brief \return Mahalanobis distance matrix for the given point index */
        inline const Eigen::Matrix4f &mahalanobis(size_t index) const {
            assert(index < mahalanobis_.size());
            return mahalanobis_[index];
        }

        /** \brief Computes rotation matrix derivative.
          * rotation matrix is obtainded from rotation angles x[3], x[4] and x[5]
          * \return d/d_rx, d/d_ry and d/d_rz respectively in g[3], g[4] and g[5]
          * param x array representing 3D transformation
          * param R rotation matrix
          * param g gradient vector
          */
        void
        computeRDerivative(const Vector6d &x, const Eigen::Matrix3d &R, Vector6d &g) const;

        /** \brief Set the rotation epsilon (maximum allowable difference between two
          * consecutive rotations) in order for an optimization to be considered as having
          * converged to the final solution.
          * \param epsilon the rotation epsilon
          */
        inline void
        setRotationEpsilon(double epsilon) { rotation_epsilon_ = epsilon; }

        /** \brief Get the rotation epsilon (maximum allowable difference between two
          * consecutive rotations) as set by the user.
          */
        inline double
        getRotationEpsilon() { return (rotation_epsilon_); }

        /** \brief Set the number of neighbors used when selecting a point neighbourhood
          * to compute covariances.
          * A higher value will bring more accurate covariance matrix but will make
          * covariances computation slower.
          * \param k the number of neighbors to use when computing covariances
          */
        void
        setCorrespondenceRandomness(int k) { k_correspondences_ = k; }

        /** \brief Get the number of neighbors used when computing covariances as set by
          * the user
          */
        int
        getCorrespondenceRandomness() { return (k_correspondences_); }

        /** set maximum number of iterations at the optimization step
          * \param[in] max maximum number of iterations for the optimizer
          */
        void
        setMaximumOptimizerIterations(int max) { max_inner_iterations_ = max; }

        ///\return maximum number of iterations at the optimization step
        int
        getMaximumOptimizerIterations() { return (max_inner_iterations_); }

    protected:

        /** \brief The number of neighbors used for covariances computation.
          * default: 20
          */
        int k_correspondences_;

        /** \brief The epsilon constant for gicp paper; this is NOT the convergence
          * tolerance
          * default: 0.001
          */
        double gicp_epsilon_;

        /** The epsilon constant for rotation error. (In GICP the transformation epsilon
          * is split in rotation part and translation part).
          * default: 2e-3
          */
        double rotation_epsilon_;

        /** \brief base transformation */
        Eigen::Matrix4f base_transformation_;

        /** \brief Temporary pointer to the source dataset. */
        const PointCloudSource *tmp_src_;

        /** \brief Temporary pointer to the target dataset. */
        const PointCloudTarget *tmp_tgt_;

        /** \brief Temporary pointer to the source dataset indices. */
        const std::vector<int> *tmp_idx_src_{};

        /** \brief Temporary pointer to the target dataset indices. */
        const std::vector<int> *tmp_idx_tgt_{};


        /** \brief Input cloud points covariances. */
        MatricesVectorPtr input_covariances_;

        /** \brief Target cloud points covariances. */
        MatricesVectorPtr target_covariances_;

        /** \brief Mahalanobis matrices holder. */
        std::vector<Eigen::Matrix4f> mahalanobis_;

        /** \brief maximum number of optimizations */
        int max_inner_iterations_;

        /** \brief compute points covariances matrices according to the K nearest
          * neighbors. K is set via setCorrespondenceRandomness() method.
          * \param cloud pointer to point cloud
          * \param tree KD tree performer for nearest neighbors search
          * \param[out] cloud_covariances covariances matrices for each point in the cloud
          */
        template<typename PointT>
        void computeCovariances(typename pcl::PointCloud<PointT>::ConstPtr cloud,
                                typename pcl::search::KdTree<PointT>::ConstPtr tree,
                                MatricesVector &cloud_covariances);

        /** \return trace of mat1^t . mat2
          * \param mat1 matrix of dimension nxm
          * \param mat2 matrix of dimension nxp
          */
        inline double
        matricesInnerProd(const Eigen::MatrixXd &mat1, const Eigen::MatrixXd &mat2) const {
            double r = 0.;
            size_t n = mat1.rows();
            // tr(mat1^t.mat2)
            for (size_t i = 0; i < n; i++)
                for (size_t j = 0; j < n; j++)
                    r += mat1(j, i) * mat2(i, j);
            return r;
        }

        /** \brief Rigid transformation computation method  with initial guess.
          * \param output the transformed input point cloud dataset using the rigid transformation found
          * \param guess the initial guess of the transformation to compute
          */
        void
        computeTransformation(PointCloudSource &output, const Eigen::Matrix4f &guess) override;

        /** \brief Search for the closest nearest neighbor of a given point.
          * \param query the point to search a nearest neighbour for
          * \param index vector of size 1 to store the index of the nearest neighbour found
          * \param distance vector of size 1 to store the distance to nearest neighbour found
          */
        inline bool
        searchForNeighbors(const PointSource &query, std::vector<int> &index, std::vector<float> &distance) const {
            int k = tree_->nearestKSearch(query, 1, index, distance);
            if (k == 0)
                return (false);
            return (true);
        }

        /// \brief compute transformation matrix from transformation matrix
        void applyState(Eigen::Matrix4f &t, const Vector6d &x) const;

        /// \brief optimization functor structure
        struct OptimizationFunctorWithIndices : public BFGSDummyFunctor<double, 6> {
            explicit OptimizationFunctorWithIndices(const GeneralizedIterativeClosestPoint *gicp)
                    : BFGSDummyFunctor<double, 6>(), gicp_(gicp) {}

            double operator()(const Vector6d &x) override;

            void df(const Vector6d &x, Vector6d &df) override;

            void fdf(const Vector6d &x, double &f, Vector6d &df) override;

            const GeneralizedIterativeClosestPoint *gicp_;
        };

        std::function<void(const pcl::PointCloud<PointSource> &cloud_src,
                           const std::vector<int> &src_indices,
                           const pcl::PointCloud<PointTarget> &cloud_tgt,
                           const std::vector<int> &tgt_indices,
                           Eigen::Matrix4f &transformation_matrix)> rigid_transformation_estimation_;
    };



////////////////////////////////////////////////////////////////////////////////////////
    template<typename PointSource, typename PointTarget>
    template<typename PointT>
    void pclomp::GeneralizedIterativeClosestPoint<PointSource, PointTarget>::computeCovariances(
            typename pcl::PointCloud<PointT>::ConstPtr cloud,
            const typename pcl::search::KdTree<PointT>::ConstPtr kdtree,
            MatricesVector &cloud_covariances) {
        if (k_correspondences_ > int(cloud->size())) {
            PCL_ERROR (
                    "[pcl::GeneralizedIterativeClosestPoint::computeCovariances] Number or points in cloud (%lu) is less than k_correspondences_ (%lu)!\n",
                    cloud->size(), k_correspondences_);
            return;
        }

        // We should never get there but who knows
        if (cloud_covariances.size() < cloud->size())
            cloud_covariances.resize(cloud->size());

        std::vector<std::vector<int>> nn_indecies_array(omp_get_max_threads());
        std::vector<std::vector<float>> nn_dist_sq_array(omp_get_max_threads());

#pragma omp parallel for
        for (int i = 0; i < cloud->size(); i++) {
            auto &nn_indecies = nn_indecies_array[omp_get_thread_num()];
            auto &nn_dist_sq = nn_dist_sq_array[omp_get_thread_num()];

            const PointT &query_point = cloud->at(i);
            Eigen::Vector3d mean = Eigen::Vector3d::Zero();
            Eigen::Matrix3d &cov = cloud_covariances[i];
            // Zero out the cov and mean
            cov.setZero();

            // Search for the K nearest neighbours
            kdtree->nearestKSearch(query_point, k_correspondences_, nn_indecies, nn_dist_sq);

            // Find the covariance matrix
            for (int j = 0; j < k_correspondences_; j++) {
                const PointT &pt = (*cloud)[nn_indecies[j]];

                mean[0] += pt.x;
                mean[1] += pt.y;
                mean[2] += pt.z;

                cov(0, 0) += pt.x * pt.x;

                cov(1, 0) += pt.y * pt.x;
                cov(1, 1) += pt.y * pt.y;

                cov(2, 0) += pt.z * pt.x;
                cov(2, 1) += pt.z * pt.y;
                cov(2, 2) += pt.z * pt.z;
            }

            mean /= static_cast<double> (k_correspondences_);
            // Get the actual covariance
            for (int k = 0; k < 3; k++)
                for (int l = 0; l <= k; l++) {
                    cov(k, l) /= static_cast<double> (k_correspondences_);
                    cov(k, l) -= mean[k] * mean[l];
                    cov(l, k) = cov(k, l);
                }

            // Compute the SVD (covariance matrix is symmetric so U = V')
            Eigen::JacobiSVD<Eigen::Matrix3d> svd(cov, Eigen::ComputeFullU);
            cov.setZero();
            Eigen::Matrix3d U = svd.matrixU();
            // Reconstitute the covariance matrix with modified singular values using the column     // vectors in V.
            for (int k = 0; k < 3; k++) {
                Eigen::Vector3d col = U.col(k);
                double v = 1.; // biggest 2 singular values replaced by 1
                if (k == 2)   // smallest singular value replaced by gicp_epsilon
                    v = gicp_epsilon_;
                cov += v * col * col.transpose();
            }
        }
    }

////////////////////////////////////////////////////////////////////////////////////////
    template<typename PointSource, typename PointTarget>
    void pclomp::GeneralizedIterativeClosestPoint<PointSource, PointTarget>::computeRDerivative(
            const Vector6d &x,
            const Eigen::Matrix3d &R,
            Vector6d &g) const {
        Eigen::Matrix3d dR_dPhi;
        Eigen::Matrix3d dR_dTheta;
        Eigen::Matrix3d dR_dPsi;

        double phi = x[3], theta = x[4], psi = x[5];

        double cphi = std::cos(phi), sphi = sin(phi);
        double ctheta = std::cos(theta), stheta = sin(theta);
        double cpsi = std::cos(psi), spsi = sin(psi);

        dR_dPhi(0, 0) = 0.;
        dR_dPhi(1, 0) = 0.;
        dR_dPhi(2, 0) = 0.;

        dR_dPhi(0, 1) = sphi * spsi + cphi * cpsi * stheta;
        dR_dPhi(1, 1) = -cpsi * sphi + cphi * spsi * stheta;
        dR_dPhi(2, 1) = cphi * ctheta;

        dR_dPhi(0, 2) = cphi * spsi - cpsi * sphi * stheta;
        dR_dPhi(1, 2) = -cphi * cpsi - sphi * spsi * stheta;
        dR_dPhi(2, 2) = -ctheta * sphi;

        dR_dTheta(0, 0) = -cpsi * stheta;
        dR_dTheta(1, 0) = -spsi * stheta;
        dR_dTheta(2, 0) = -ctheta;

        dR_dTheta(0, 1) = cpsi * ctheta * sphi;
        dR_dTheta(1, 1) = ctheta * sphi * spsi;
        dR_dTheta(2, 1) = -sphi * stheta;

        dR_dTheta(0, 2) = cphi * cpsi * ctheta;
        dR_dTheta(1, 2) = cphi * ctheta * spsi;
        dR_dTheta(2, 2) = -cphi * stheta;

        dR_dPsi(0, 0) = -ctheta * spsi;
        dR_dPsi(1, 0) = cpsi * ctheta;
        dR_dPsi(2, 0) = 0.;

        dR_dPsi(0, 1) = -cphi * cpsi - sphi * spsi * stheta;
        dR_dPsi(1, 1) = -cphi * spsi + cpsi * sphi * stheta;
        dR_dPsi(2, 1) = 0.;

        dR_dPsi(0, 2) = cpsi * sphi - cphi * spsi * stheta;
        dR_dPsi(1, 2) = sphi * spsi + cphi * cpsi * stheta;
        dR_dPsi(2, 2) = 0.;

        g[3] = matricesInnerProd(dR_dPhi, R);
        g[4] = matricesInnerProd(dR_dTheta, R);
        g[5] = matricesInnerProd(dR_dPsi, R);
    }

////////////////////////////////////////////////////////////////////////////////////////
    template<typename PointSource, typename PointTarget>
    void pclomp::GeneralizedIterativeClosestPoint<PointSource, PointTarget>::estimateRigidTransformationBFGS(
            const PointCloudSource &cloud_src,
            const std::vector<int> &indices_src,
            const PointCloudTarget &cloud_tgt,
            const std::vector<int> &indices_tgt,
            Eigen::Matrix4f &transformation_matrix) {
        if (indices_src.size() < 4)     // need at least 4 samples
        {
            PCL_THROW_EXCEPTION (pcl::NotEnoughPointsException,
                                 "[pcl::GeneralizedIterativeClosestPoint::estimateRigidTransformationBFGS] Need at least 4 points to estimate a transform! Source and target have "
                                         << indices_src.size() << " points!");
            return;
        }
        // Set the initial solution
        Vector6d x = Vector6d::Zero();
        x[0] = transformation_matrix(0, 3);
        x[1] = transformation_matrix(1, 3);
        x[2] = transformation_matrix(2, 3);
        x[3] = std::atan2(transformation_matrix(2, 1), transformation_matrix(2, 2));
        x[4] = asin(-transformation_matrix(2, 0));
        x[5] = std::atan2(transformation_matrix(1, 0), transformation_matrix(0, 0));

        // Set temporary pointers
        tmp_src_ = &cloud_src;
        tmp_tgt_ = &cloud_tgt;
        tmp_idx_src_ = &indices_src;
        tmp_idx_tgt_ = &indices_tgt;

        // Optimize using forward-difference approximation LM
        const double gradient_tol = 1e-2;
        OptimizationFunctorWithIndices functor(this);
        BFGS<OptimizationFunctorWithIndices> bfgs(functor);
        bfgs.parameters.sigma = 0.01;
        bfgs.parameters.rho = 0.01;
        bfgs.parameters.tau1 = 9;
        bfgs.parameters.tau2 = 0.05;
        bfgs.parameters.tau3 = 0.5;
        bfgs.parameters.order = 3;

        int inner_iterations_ = 0;
        int result = bfgs.minimizeInit(x);
        result = BFGSSpace::Running;
        do {
            inner_iterations_++;
            result = bfgs.minimizeOneStep(x);
            if (result) {
                break;
            }
            result = bfgs.testGradient(gradient_tol);
        } while (result == BFGSSpace::Running && inner_iterations_ < max_inner_iterations_);
        if (result == BFGSSpace::NoProgress || result == BFGSSpace::Success || inner_iterations_ == max_inner_iterations_) {
            PCL_DEBUG ("[pcl::registration::TransformationEstimationBFGS::estimateRigidTransformation]");
            PCL_DEBUG ("BFGS solver finished with exit code %i \n", result);
            transformation_matrix.setIdentity();
            applyState(transformation_matrix, x);
        } else PCL_THROW_EXCEPTION(pcl::SolverDidntConvergeException,
                                   "[pcl::" << getClassName()
                                            << "::TransformationEstimationBFGS::estimateRigidTransformation] BFGS solver didn't converge!");
    }

////////////////////////////////////////////////////////////////////////////////////////
    template<typename PointSource, typename PointTarget>
    inline double
    pclomp::GeneralizedIterativeClosestPoint<PointSource, PointTarget>::OptimizationFunctorWithIndices::operator()(
            const Vector6d &x) {
        Eigen::Matrix4f transformation_matrix = gicp_->base_transformation_;
        gicp_->applyState(transformation_matrix, x);
        double f = 0;
        std::vector<double> f_array(omp_get_max_threads(), 0.0);

        int m = static_cast<int> (gicp_->tmp_idx_src_->size());
#pragma omp parallel for
        for (int i = 0; i < m; ++i) {
            // The last coordinate, p_src[3] is guaranteed to be set to 1.0 in registration.hpp
            pcl::Vector4fMapConst p_src = gicp_->tmp_src_->points[(*gicp_->tmp_idx_src_)[i]].getVector4fMap();
            // The last coordinate, p_tgt[3] is guaranteed to be set to 1.0 in registration.hpp
            pcl::Vector4fMapConst p_tgt = gicp_->tmp_tgt_->points[(*gicp_->tmp_idx_tgt_)[i]].getVector4fMap();
            // Estimate the distance (cost function)
            // The last coordiante is still guaranteed to be set to 1.0
            // Eigen::AlignedVector3<double> res(pp[0] - p_tgt[0], pp[1] - p_tgt[1], pp[2] - p_tgt[2]);
            // Eigen::AlignedVector3<double> temp(gicp_->mahalanobis((*gicp_->tmp_idx_src_)[i]) * res);
            Eigen::Vector4f res = transformation_matrix * p_src - p_tgt;
            Eigen::Matrix4f maha = gicp_->mahalanobis((*gicp_->tmp_idx_src_)[i]);
            // Eigen::Vector4d temp(maha * res);
            // increment= res'*temp/num_matches = temp'*M*temp/num_matches (we postpone 1/num_matches after the loop closes)
            // double ret = double(res.transpose() * temp);
            double ret = res.dot(maha * res);
            f_array[omp_get_thread_num()] += ret;
        }
        f = std::accumulate(f_array.begin(), f_array.end(), 0.0);
        return f / m;
    }

////////////////////////////////////////////////////////////////////////////////////////
    template<typename PointSource, typename PointTarget>
    inline void pclomp::GeneralizedIterativeClosestPoint<PointSource, PointTarget>::OptimizationFunctorWithIndices::df(
            const Vector6d &x, Vector6d &g) {
        Eigen::Matrix4f transformation_matrix = gicp_->base_transformation_;
        gicp_->applyState(transformation_matrix, x);
        //Eigen::Vector3d g_t = g.head<3> ();
        std::vector<Eigen::Matrix4d, Eigen::aligned_allocator<Eigen::Matrix4d>> R_array(omp_get_max_threads());
        std::vector<Eigen::Vector4d, Eigen::aligned_allocator<Eigen::Vector4d>> g_array(omp_get_max_threads());

        for (int i = 0; i < R_array.size(); i++) {
            R_array[i].setZero();
            g_array[i].setZero();
        }

        int m = static_cast<int> (gicp_->tmp_idx_src_->size());
#pragma omp parallel for
        for (int i = 0; i < m; ++i) {
            // The last coordinate, p_src[3] is guaranteed to be set to 1.0 in registration.hpp
            pcl::Vector4fMapConst p_src = gicp_->tmp_src_->points[(*gicp_->tmp_idx_src_)[i]].getVector4fMap();
            // The last coordinate, p_tgt[3] is guaranteed to be set to 1.0 in registration.hpp
            pcl::Vector4fMapConst p_tgt = gicp_->tmp_tgt_->points[(*gicp_->tmp_idx_tgt_)[i]].getVector4fMap();

            Eigen::Vector4f pp(transformation_matrix * p_src);
            // The last coordiante is still guaranteed to be set to 1.0
            Eigen::Vector4d res(pp[0] - p_tgt[0], pp[1] - p_tgt[1], pp[2] - p_tgt[2], 0.0);
            // temp = M*res

            Eigen::Matrix4d maha = gicp_->mahalanobis((*gicp_->tmp_idx_src_)[i]).template cast<double>();

            Eigen::Vector4d temp(maha * res);
            // Increment translation gradient
            // g.head<3> ()+= 2*M*res/num_matches (we postpone 2/num_matches after the loop closes)
            // Increment rotation gradient
            pp = gicp_->base_transformation_ * p_src;

            Eigen::Vector4d p_src3(pp[0], pp[1], pp[2], 0.0);
            g_array[omp_get_thread_num()] += temp;
            R_array[omp_get_thread_num()] += p_src3 * temp.transpose();
        }

        g.setZero();
        Eigen::Matrix4d R = Eigen::Matrix4d::Zero();
        for (int i = 0; i < R_array.size(); i++) {
            R += R_array[i];
            g.head<3>() += g_array[i].head<3>();
        }

        g.head<3>() *= 2.0 / m;
        R *= 2.0 / m;

        gicp_->computeRDerivative(x, R.block<3, 3>(0, 0), g);
    }

////////////////////////////////////////////////////////////////////////////////////////
    template<typename PointSource, typename PointTarget>
    inline void pclomp::GeneralizedIterativeClosestPoint<PointSource, PointTarget>::OptimizationFunctorWithIndices::fdf(
            const Vector6d &x, double &f, Vector6d &g) {
        Eigen::Matrix4f transformation_matrix = gicp_->base_transformation_;
        gicp_->applyState(transformation_matrix, x);
        f = 0;
        g.setZero();
        Eigen::Matrix3d R = Eigen::Matrix3d::Zero();
        const int m = static_cast<const int> (gicp_->tmp_idx_src_->size());
        for (int i = 0; i < m; ++i) {
            // The last coordinate, p_src[3] is guaranteed to be set to 1.0 in registration.hpp
            pcl::Vector4fMapConst p_src = gicp_->tmp_src_->points[(*gicp_->tmp_idx_src_)[i]].getVector4fMap();
            // The last coordinate, p_tgt[3] is guaranteed to be set to 1.0 in registration.hpp
            pcl::Vector4fMapConst p_tgt = gicp_->tmp_tgt_->points[(*gicp_->tmp_idx_tgt_)[i]].getVector4fMap();
            Eigen::Vector4f pp(transformation_matrix * p_src);
            // The last coordinate is still guaranteed to be set to 1.0
            Eigen::Vector3d res(pp[0] - p_tgt[0], pp[1] - p_tgt[1], pp[2] - p_tgt[2]);
            // temp = M*res
            Eigen::Vector3d temp(gicp_->mahalanobis((*gicp_->tmp_idx_src_)[i]).template block<3, 3>(0, 0).
                    template cast<double>() * res);
            // Increment total error
            f += double(res.transpose() * temp);
            // Increment translation gradient
            // g.head<3> ()+= 2*M*res/num_matches (we postpone 2/num_matches after the loop closes)
            g.head<3>() += temp;
            pp = gicp_->base_transformation_ * p_src;
            Eigen::Vector3d p_src3(pp[0], pp[1], pp[2]);
            // Increment rotation gradient
            R += p_src3 * temp.transpose();
        }
        f /= double(m);
        g.head<3>() *= double(2.0 / m);
        R *= 2.0 / m;
        gicp_->computeRDerivative(x, R, g);
    }

////////////////////////////////////////////////////////////////////////////////////////
    template<typename PointSource, typename PointTarget>
    inline void pclomp::GeneralizedIterativeClosestPoint<PointSource, PointTarget>::computeTransformation(
            PointCloudSource &output,
            const Eigen::Matrix4f &guess) {
        pcl::IterativeClosestPoint<PointSource, PointTarget>::initComputeReciprocal();
        using namespace std;
        // Difference between consecutive transforms
        double delta = 0;
        // Get the size of the target
        const size_t N = indices_->size();
        // Set the mahalanobis matrices to identity
        mahalanobis_.resize(N, Eigen::Matrix4f::Identity());
        // Compute target cloud covariance matrices
        if ((!target_covariances_) || (target_covariances_->empty())) {
            target_covariances_.reset(new MatricesVector);
            computeCovariances<PointTarget>(target_, tree_, *target_covariances_);
        }
        // Compute input cloud covariance matrices
        if ((!input_covariances_) || (input_covariances_->empty())) {
            input_covariances_.reset(new MatricesVector);
            computeCovariances<PointSource>(input_, tree_reciprocal_, *input_covariances_);
        }

        base_transformation_ = Eigen::Matrix4f::Identity();
        nr_iterations_ = 0;
        converged_ = false;
        double dist_threshold = corr_dist_threshold_ * corr_dist_threshold_;
        pcl::transformPointCloud(output, output, guess);

        std::vector<std::vector<int>> nn_indices_array(omp_get_max_threads());
        std::vector<std::vector<float>> nn_dists_array(omp_get_max_threads());
        for (auto &nn_indices: nn_indices_array) { nn_indices.resize(1); }
        for (auto &nn_dists: nn_dists_array) { nn_dists.resize(1); }

        while (!converged_) {
            std::atomic<size_t> cnt;
            cnt = 0;
            std::vector<int> source_indices(indices_->size());
            std::vector<int> target_indices(indices_->size());

            // guess corresponds to base_t and transformation_ to t
            Eigen::Matrix4d transform_R = Eigen::Matrix4d::Zero();
            for (size_t i = 0; i < 4; i++)
                for (size_t j = 0; j < 4; j++)
                    for (size_t k = 0; k < 4; k++)
                        transform_R(i, j) += double(transformation_(i, k)) * double(guess(k, j));

            const Eigen::Matrix3d R = transform_R.topLeftCorner<3, 3>();

#pragma omp parallel for
            for (int i = 0; i < N; i++) {
                auto &nn_indices = nn_indices_array[omp_get_thread_num()];
                auto &nn_dists = nn_dists_array[omp_get_thread_num()];

                PointSource query = output[i];
                query.getVector4fMap() = transformation_ * query.getVector4fMap();

                if (!searchForNeighbors(query, nn_indices, nn_dists)) {
                    PCL_ERROR (
                            "[pcl::%s::computeTransformation] Unable to find a nearest neighbor in the target dataset for point %d in the source!\n",
                            getClassName().c_str(), (*indices_)[i]);
                    continue;
                }

                // Check if the distance to the nearest neighbor is smaller than the user imposed threshold
                if (nn_dists[0] < dist_threshold) {
                    const Eigen::Matrix3d &C1 = (*input_covariances_)[i];
                    const Eigen::Matrix3d &C2 = (*target_covariances_)[nn_indices[0]];
                    Eigen::Matrix4f &M_ = mahalanobis_[i];
                    M_.setZero();

                    Eigen::Matrix3d M = M_.block<3, 3>(0, 0).cast<double>();
                    // M = R*C1
                    M = R * C1;
                    // temp = M*R' + C2 = R*C1*R' + C2
                    Eigen::Matrix3d temp = M * R.transpose();
                    temp += C2;
                    // M = temp^-1
                    M = temp.inverse();
                    M_.block<3, 3>(0, 0) = M.cast<float>();
                    int c = cnt++;
                    source_indices[c] = static_cast<int> (i);
                    target_indices[c] = nn_indices[0];
                }
            }
            // Resize to the actual number of valid correspondences
            source_indices.resize(cnt);
            target_indices.resize(cnt);

            std::vector<std::pair<int, int>> indices(source_indices.size());
            for (int i = 0; i < source_indices.size(); i++) {
                indices[i].first = source_indices[i];
                indices[i].second = target_indices[i];
            }

            std::sort(indices.begin(), indices.end(), [=](const std::pair<int, int> &lhs, const std::pair<int, int> &rhs) {
                return lhs.first < rhs.first;
            });

            for (int i = 0; i < source_indices.size(); i++) {
                source_indices[i] = indices[i].first;
                target_indices[i] = indices[i].second;
            }

            /* optimize transformation using the current assignment and Mahalanobis metrics*/
            previous_transformation_ = transformation_;
            //optimization right here
            try {
                rigid_transformation_estimation_(output, source_indices, *target_, target_indices, transformation_);
                /* compute the delta from this iteration */
                delta = 0.;
                for (int k = 0; k < 4; k++) {
                    for (int l = 0; l < 4; l++) {
                        double ratio = 1;
                        if (k < 3 && l < 3) // rotation part of the transform
                            ratio = 1. / rotation_epsilon_;
                        else
                            ratio = 1. / transformation_epsilon_;
                        double c_delta = ratio * std::abs(previous_transformation_(k, l) - transformation_(k, l));
                        if (c_delta > delta)
                            delta = c_delta;
                    }
                }
            }
            catch (pcl::PCLException &e) {
                PCL_DEBUG ("[pcl::%s::computeTransformation] Optimization issue %s\n", getClassName().c_str(), e.what());
                break;
            }
            nr_iterations_++;
            // Check for convergence
            if (nr_iterations_ >= max_iterations_ || delta < 1) {
                converged_ = true;
                previous_transformation_ = transformation_;
                PCL_DEBUG (
                        "[pcl::%s::computeTransformation] Convergence reached. Number of iterations: %d out of %d. Transformation difference: %f\n",
                        getClassName().c_str(), nr_iterations_, max_iterations_,
                        (transformation_ - previous_transformation_).array().abs().sum());
            } else
                PCL_DEBUG ("[pcl::%s::computeTransformation] Convergence failed\n", getClassName().c_str());
        }
        final_transformation_ = previous_transformation_ * guess;

        // Transform the point cloud
        pcl::transformPointCloud(*input_, output, final_transformation_);
    }

    template<typename PointSource, typename PointTarget>
    void pclomp::GeneralizedIterativeClosestPoint<PointSource, PointTarget>::applyState(
            Eigen::Matrix4f &t,
            const Vector6d &x) const {
        // !!! CAUTION Stanford GICP uses the Z Y X euler angles convention
        Eigen::Matrix3f R;
        R = Eigen::AngleAxisf(static_cast<float> (x[5]), Eigen::Vector3f::UnitZ())
            * Eigen::AngleAxisf(static_cast<float> (x[4]), Eigen::Vector3f::UnitY())
            * Eigen::AngleAxisf(static_cast<float> (x[3]), Eigen::Vector3f::UnitX());
        t.topLeftCorner<3, 3>().matrix() = R * t.topLeftCorner<3, 3>().matrix();
        Eigen::Vector4f T(static_cast<float> (x[0]), static_cast<float> (x[1]), static_cast<float> (x[2]), 0.0f);
        t.col(3) += T;
    }

}

#endif  //PCL_GICP_OMP_H_
