#ifndef PCL_VOXEL_GRID_COVARIANCE_OMP_H_
#define PCL_VOXEL_GRID_COVARIANCE_OMP_H_

#include "Eigen/Cholesky"
#include "Eigen/Dense"
#include "map"
#include "pcl/common/common.h"
#include "pcl/filters/boost.h"
#include "pcl/filters/voxel_grid.h"
#include "pcl/kdtree/kdtree_flann.h"
#include "pcl/point_types.h"
#include "pclomp/voxel_grid_covariance_omp.hpp"
#include "unordered_map"

namespace pclomp {
/** \brief A searchable voxel strucure containing the mean and covariance of the data.
 * \note For more information please see
 * <b>Magnusson, M. (2009). The Three-Dimensional Normal-Distributions Transform —
 * an Efﬁcient Representation for Registration, Surface Analysis, and Loop Detection.
 * PhD thesis, Orebro University. Orebro Studies in Technology 36</b>
 * \author Brian Okorn (Space and Naval Warfare Systems Center Pacific)
 */
template <typename PointT>
class VoxelGridCovariance : public pcl::VoxelGrid<PointT> {
protected:
    using pcl::VoxelGrid<PointT>::filter_name_;
    using pcl::VoxelGrid<PointT>::getClassName;
    using pcl::VoxelGrid<PointT>::input_;
    using pcl::VoxelGrid<PointT>::indices_;
    using pcl::VoxelGrid<PointT>::filter_limit_negative_;
    using pcl::VoxelGrid<PointT>::filter_limit_min_;
    using pcl::VoxelGrid<PointT>::filter_limit_max_;
    using pcl::VoxelGrid<PointT>::filter_field_name_;

    using pcl::VoxelGrid<PointT>::downsample_all_data_;
    using pcl::VoxelGrid<PointT>::leaf_layout_;
    using pcl::VoxelGrid<PointT>::save_leaf_layout_;
    using pcl::VoxelGrid<PointT>::leaf_size_;
    using pcl::VoxelGrid<PointT>::min_b_;
    using pcl::VoxelGrid<PointT>::max_b_;
    using pcl::VoxelGrid<PointT>::inverse_leaf_size_;
    using pcl::VoxelGrid<PointT>::div_b_;
    using pcl::VoxelGrid<PointT>::divb_mul_;

    typedef typename pcl::traits::fieldList<PointT>::type FieldList;
    typedef typename pcl::Filter<PointT>::PointCloud PointCloud;
    typedef typename PointCloud::Ptr PointCloudPtr;
    typedef typename PointCloud::ConstPtr PointCloudConstPtr;

public:
    typedef boost::shared_ptr<pcl::VoxelGrid<PointT>> Ptr;
    typedef boost::shared_ptr<const pcl::VoxelGrid<PointT>> ConstPtr;

    /** \brief Simple structure to hold a centroid, covarince and the number of points in a leaf.
     * Inverse covariance, eigen vectors and engen values are precomputed. */
    struct Leaf {
        /** \brief Constructor.
         * Sets \ref nr_points, \ref icov_, \ref mean_ and \ref evals_ to 0 and \ref cov_ and \ref
         * evecs_ to the identity matrix
         */
        Leaf()
            : nr_points(0),
              mean_(Eigen::Vector3d::Zero()),
              centroid(),
              cov_(Eigen::Matrix3d::Identity()),
              icov_(Eigen::Matrix3d::Zero()),
              evecs_(Eigen::Matrix3d::Identity()),
              evals_(Eigen::Vector3d::Zero()) {}

        /** \brief Get the voxel covariance.
         * \return covariance matrix
         */
        Eigen::Matrix3d getCov() const { return (cov_); }

        /** \brief Get the inverse of the voxel covariance.
         * \return inverse covariance matrix
         */
        Eigen::Matrix3d getInverseCov() const { return (icov_); }

        /** \brief Get the voxel centroid.
         * \return centroid
         */
        Eigen::Vector3d getMean() const { return (mean_); }

        /** \brief Get the eigen vectors of the voxel covariance.
         * \note Order corresponds with \ref getEvals
         * \return matrix whose columns contain eigen vectors
         */
        Eigen::Matrix3d getEvecs() const { return (evecs_); }

        /** \brief Get the eigen values of the voxel covariance.
         * \note Order corresponds with \ref getEvecs
         * \return vector of eigen values
         */
        Eigen::Vector3d getEvals() const { return (evals_); }

        /** \brief Get the number of points contained by this voxel.
         * \return number of points
         */
        int getPointCount() const { return (nr_points); }

        /** \brief Number of points contained by voxel */
        int nr_points;

        /** \brief 3D voxel centroid */
        Eigen::Vector3d mean_;

        /** \brief Nd voxel centroid
         * \note Differs from \ref mean_ when color data is used
         */
        Eigen::VectorXf centroid;

        /** \brief Voxel covariance matrix */
        Eigen::Matrix3d cov_;

        /** \brief Inverse of voxel covariance matrix */
        Eigen::Matrix3d icov_;

        /** \brief Eigen vectors of voxel covariance matrix */
        Eigen::Matrix3d evecs_;

        /** \brief Eigen values of voxel covariance matrix */
        Eigen::Vector3d evals_;

        /** \brief Points inside the cell */
        pcl::PointCloud<PointT> pointList_;
    };

    /** \brief Pointer to VoxelGridCovariance leaf structure */
    typedef Leaf *LeafPtr;

    /** \brief Const pointer to VoxelGridCovariance leaf structure */
    typedef const Leaf *LeafConstPtr;

    typedef std::map<size_t, Leaf> Map;

public:
    /** \brief Constructor.
     * Sets \ref leaf_size_ to 0 and \ref searchable_ to false.
     */
    VoxelGridCovariance()
        : searchable_(true),
          min_points_per_voxel_(6),
          min_covar_eigvalue_mult_(0.01),
          leaves_(),
          voxel_centroids_(),
          voxel_centroids_leaf_indices_(),
          kdtree_() {
        downsample_all_data_ = false;
        save_leaf_layout_ = false;
        leaf_size_.setZero();
        min_b_.setZero();
        max_b_.setZero();
        filter_name_ = "VoxelGridCovariance";
    }

    /** \brief Set the minimum number of points required for a cell to be used (must be 3 or greater
     * for covariance calculation).
     * \param[in] min_points_per_voxel the minimum number of points for required for a voxel to be
     * used
     */
    inline void setMinPointPerVoxel(int min_points_per_voxel) {
        if (min_points_per_voxel > 2) {
            min_points_per_voxel_ = min_points_per_voxel;
        } else {
            PCL_WARN(
                "%s: Covariance calculation requires at least 3 points, setting Min Point per "
                "Voxel to 3 ",
                this->getClassName().c_str());
            min_points_per_voxel_ = 3;
        }
    }

    /** \brief Get the minimum number of points required for a cell to be used.
     * \return the minimum number of points for required for a voxel to be used
     */
    inline int getMinPointPerVoxel() { return min_points_per_voxel_; }

    /** \brief Set the minimum allowable ratio between eigenvalues to prevent singular covariance
     * matrices.
     * \param[in] min_covar_eigvalue_mult the minimum allowable ratio between eigenvalues
     */
    inline void setCovEigValueInflationRatio(double min_covar_eigvalue_mult) {
        min_covar_eigvalue_mult_ = min_covar_eigvalue_mult;
    }

    /** \brief Get the minimum allowable ratio between eigenvalues to prevent singular covariance
     * matrices.
     * \return the minimum allowable ratio between eigenvalues
     */
    inline double getCovEigValueInflationRatio() { return min_covar_eigvalue_mult_; }

    /** \brief Filter cloud and initializes voxel structure.
     * \param[out] output cloud containing centroids of voxels containing a sufficient number of
     * points
     * \param[in] searchable flag if voxel structure is searchable, if true then kdtree is built
     */
    inline void filter(PointCloud &output, bool searchable = false) {
        searchable_ = searchable;
        applyFilter(output);

        voxel_centroids_ = PointCloudPtr(new PointCloud(output));

        if (searchable_ && voxel_centroids_->size() > 0) {
            // Initiates kdtree of the centroids of voxels containing a sufficient number of points
            kdtree_.setInputCloud(voxel_centroids_);
        }
    }

    /** \brief Initializes voxel structure.
     * \param[in] searchable flag if voxel structure is searchable, if true then kdtree is built
     */
    inline void filter(bool searchable = false) {
        searchable_ = searchable;
        voxel_centroids_ = PointCloudPtr(new PointCloud);
        applyFilter(*voxel_centroids_);

        if (searchable_ && voxel_centroids_->size() > 0) {
            // Initiates kdtree of the centroids of voxels containing a sufficient number of points
            kdtree_.setInputCloud(voxel_centroids_);
        }
    }

    /** \brief Get the voxel containing point p.
     * \param[in] index the index of the leaf structure node
     * \return const pointer to leaf structure
     */
    inline LeafConstPtr getLeaf(int index) {
        auto leaf_iter = leaves_.find(index);
        if (leaf_iter != leaves_.end()) {
            LeafConstPtr ret(&(leaf_iter->second));
            return ret;
        } else
            return NULL;
    }

    /** \brief Get the voxel containing point p.
     * \param[in] p the point to get the leaf structure at
     * \return const pointer to leaf structure
     */
    inline LeafConstPtr getLeaf(PointT &p) {
        // Generate index associated with p
        int ijk0 = static_cast<int>(floor(p.x * inverse_leaf_size_[0]) - min_b_[0]);
        int ijk1 = static_cast<int>(floor(p.y * inverse_leaf_size_[1]) - min_b_[1]);
        int ijk2 = static_cast<int>(floor(p.z * inverse_leaf_size_[2]) - min_b_[2]);

        // Compute the centroid leaf index
        int idx = ijk0 * divb_mul_[0] + ijk1 * divb_mul_[1] + ijk2 * divb_mul_[2];

        // Find leaf associated with index
        auto leaf_iter = leaves_.find(idx);
        if (leaf_iter != leaves_.end()) {
            // If such a leaf exists return the pointer to the leaf structure
            LeafConstPtr ret(&(leaf_iter->second));
            return ret;
        } else
            return NULL;
    }

    /** \brief Get the voxel containing point p.
     * \param[in] p the point to get the leaf structure at
     * \return const pointer to leaf structure
     */
    inline LeafConstPtr getLeaf(Eigen::Vector3f &p) {
        // Generate index associated with p
        int ijk0 = static_cast<int>(floor(p[0] * inverse_leaf_size_[0]) - min_b_[0]);
        int ijk1 = static_cast<int>(floor(p[1] * inverse_leaf_size_[1]) - min_b_[1]);
        int ijk2 = static_cast<int>(floor(p[2] * inverse_leaf_size_[2]) - min_b_[2]);

        // Compute the centroid leaf index
        int idx = ijk0 * divb_mul_[0] + ijk1 * divb_mul_[1] + ijk2 * divb_mul_[2];

        // Find leaf associated with index
        auto leaf_iter = leaves_.find(idx);
        if (leaf_iter != leaves_.end()) {
            // If such a leaf exists return the pointer to the leaf structure
            LeafConstPtr ret(&(leaf_iter->second));
            return ret;
        } else
            return NULL;
    }

    /** \brief Get the voxels surrounding point p, not including the voxel contating point p.
     * \note Only voxels containing a sufficient number of points are used (slower than radius
     * search in practice).
     * \param[in] reference_point the point to get the leaf structure at
     * \param[out] neighbors
     * \return number of neighbors found
     */
    int getNeighborhoodAtPoint(const Eigen::MatrixXi &,
                               const PointT &reference_point,
                               std::vector<LeafConstPtr> &neighbors) const;

    int getNeighborhoodAtPoint(const PointT &reference_point,
                               std::vector<LeafConstPtr> &neighbors) const;

    int getNeighborhoodAtPoint7(const PointT &reference_point,
                                std::vector<LeafConstPtr> &neighbors) const;

    int getNeighborhoodAtPoint1(const PointT &reference_point,
                                std::vector<LeafConstPtr> &neighbors) const;

    /** \brief Get the leaf structure map
     * \return a map contataining all leaves
     */
    inline const Map &getLeaves() const { return leaves_; }

    /** \brief Get a pointcloud containing the voxel centroids
     * \note Only voxels containing a sufficient number of points are used.
     * \return a map contataining all leaves
     */
    inline PointCloudPtr getCentroids() const { return voxel_centroids_; }

    /** \brief Get a cloud to visualize each voxels normal distribution.
     * \param[out] cell_cloud a cloud created by sampling the normal distributions of each voxel
     */
    void getDisplayCloud(pcl::PointCloud<pcl::PointXYZ> &cell_cloud);

    /** \brief Search for the k-nearest occupied voxels for the given query point.
     * \note Only voxels containing a sufficient number of points are used.
     * \param[in] point the given query point
     * \param[in] k the number of neighbors to search for
     * \param[out] k_leaves the resultant leaves of the neighboring points
     * \param[out] k_sqr_distances the resultant squared distances to the neighboring points
     * \return number of neighbors found
     */
    int nearestKSearch(const PointT &point,
                       int k,
                       std::vector<LeafConstPtr> &k_leaves,
                       std::vector<float> &k_sqr_distances) {
        k_leaves.clear();

        // Check if kdtree has been built
        if (!searchable_) {
            PCL_WARN("%s: Not Searchable", this->getClassName().c_str());
            return 0;
        }

        // Find k-nearest neighbors in the occupied voxel centroid cloud
        std::vector<int> k_indices;
        k = kdtree_.nearestKSearch(point, k, k_indices, k_sqr_distances);

        // Find leaves corresponding to neighbors
        k_leaves.reserve(k);
        for (int &k_indice : k_indices) {
            k_leaves.push_back(&leaves_[voxel_centroids_leaf_indices_[k_indice]]);
        }
        return k;
    }

    /** \brief Search for the k-nearest occupied voxels for the given query point.
     * \note Only voxels containing a sufficient number of points are used.
     * \param[in] cloud the given query point
     * \param[in] index the index
     * \param[in] k the number of neighbors to search for
     * \param[out] k_leaves the resultant leaves of the neighboring points
     * \param[out] k_sqr_distances the resultant squared distances to the neighboring points
     * \return number of neighbors found
     */
    inline int nearestKSearch(const PointCloud &cloud,
                              int index,
                              int k,
                              std::vector<LeafConstPtr> &k_leaves,
                              std::vector<float> &k_sqr_distances) {
        if (index >= static_cast<int>(cloud.points.size()) || index < 0) return (0);
        return (nearestKSearch(cloud.points[index], k, k_leaves, k_sqr_distances));
    }

    /** \brief Search for all the nearest occupied voxels of the query point in a given radius.
     * \note Only voxels containing a sufficient number of points are used.
     * \param[in] point the given query point
     * \param[in] radius the radius of the sphere bounding all of p_q's neighbors
     * \param[out] k_leaves the resultant leaves of the neighboring points
     * \param[out] k_sqr_distances the resultant squared distances to the neighboring points
     * \param[in] max_nn
     * \return number of neighbors found
     */
    int radiusSearch(const PointT &point,
                     double radius,
                     std::vector<LeafConstPtr> &k_leaves,
                     std::vector<float> &k_sqr_distances,
                     unsigned int max_nn = 0) const {
        k_leaves.clear();

        // Check if kdtree has been built
        if (!searchable_) {
            PCL_WARN("%s: Not Searchable", this->getClassName().c_str());
            return 0;
        }

        // Find neighbors within radius in the occupied voxel centroid cloud
        std::vector<int> k_indices;
        int k = kdtree_.radiusSearch(point, radius, k_indices, k_sqr_distances, max_nn);

        // Find leaves corresponding to neighbors
        k_leaves.reserve(k);
        for (int &k_indice : k_indices) {
            auto leaf = leaves_.find(voxel_centroids_leaf_indices_[k_indice]);
            if (leaf == leaves_.end()) {
                std::cerr << "error : could not find the leaf corresponding to the voxel"
                          << std::endl;
                std::cin.ignore(1);
            }
            k_leaves.push_back(&(leaf->second));
        }
        return k;
    }

    int radiusSearch(const PointT &point,
                     double radius,
                     std::vector<LeafConstPtr> &k_leaves,
                     std::vector<int> &k_indices,
                     std::vector<float> &k_sqr_distances,
                     unsigned int max_nn = 0) const {
        k_leaves.clear();

        // Check if kdtree has been built
        if (!searchable_) {
            PCL_WARN("%s: Not Searchable", this->getClassName().c_str());
            return 0;
        }

        // Find neighbors within radius in the occupied voxel centroid cloud
        int k = kdtree_.radiusSearch(point, radius, k_indices, k_sqr_distances, max_nn);

        // Find leaves corresponding to neighbors
        k_leaves.reserve(k);
        for (int &k_indice : k_indices) {
            auto leaf = leaves_.find(voxel_centroids_leaf_indices_[k_indice]);
            if (leaf == leaves_.end()) {
                std::cerr << "error : could not find the leaf corresponding to the voxel"
                          << std::endl;
                std::cin.ignore(1);
            }
            k_leaves.push_back(&(leaf->second));
        }
        return k;
    }

    /** \brief Search for all the nearest occupied voxels of the query point in a given radius.
     * \note Only voxels containing a sufficient number of points are used.
     * \param[in] cloud the given query point
     * \param[in] index a valid index in cloud representing a valid (i.e., finite) query point
     * \param[in] radius the radius of the sphere bounding all of p_q's neighbors
     * \param[out] k_leaves the resultant leaves of the neighboring points
     * \param[out] k_sqr_distances the resultant squared distances to the neighboring points
     * \param[in] max_nn
     * \return number of neighbors found
     */
    inline int radiusSearch(const PointCloud &cloud,
                            int index,
                            double radius,
                            std::vector<LeafConstPtr> &k_leaves,
                            std::vector<float> &k_sqr_distances,
                            unsigned int max_nn = 0) const {
        if (index >= static_cast<int>(cloud.points.size()) || index < 0) return (0);
        return (radiusSearch(cloud.points[index], radius, k_leaves, k_sqr_distances, max_nn));
    }

protected:
    /** \brief Filter cloud and initializes voxel structure.
     * \param[out] output cloud containing centroids of voxels containing a sufficient number of
     * points
     */
    void applyFilter(PointCloud &output);

    /** \brief Flag to determine if voxel structure is searchable. */
    bool searchable_;

    /** \brief Minimum points contained with in a voxel to allow it to be useable. */
    int min_points_per_voxel_;

    /** \brief Minimum allowable ratio between eigenvalues to prevent singular covariance matrices.
     */
    double min_covar_eigvalue_mult_;

    /** \brief Voxel structure containing all leaf nodes (includes voxels with less than a
     * sufficient number of points). */
    Map leaves_;

    /** \brief Point cloud containing centroids of voxels containing atleast minimum number of
     * points. */
    PointCloudPtr voxel_centroids_;

    /** \brief Indices of leaf structurs associated with each point in \ref voxel_centroids_ (used
     * for searching). */
    std::vector<int> voxel_centroids_leaf_indices_;

    /** \brief KdTree generated using \ref voxel_centroids_ (used for searching). */
    pcl::KdTreeFLANN<PointT> kdtree_;
};

//////////////////////////////////////////////////////////////////////////////////////////
template <typename PointT>
void pclomp::VoxelGridCovariance<PointT>::applyFilter(PointCloud &output) {
    voxel_centroids_leaf_indices_.clear();

    // Has the input dataset been set already?
    if (!input_) {
        PCL_WARN("[pcl::%s::applyFilter] No input dataset given!\n", getClassName().c_str());
        output.width = output.height = 0;
        output.points.clear();
        return;
    }

    // Copy the header (and thus the frame_id) + allocate enough space for points
    output.height = 1;       // down sampling breaks the organized structure
    output.is_dense = true;  // we filter out invalid points
    output.points.clear();

    Eigen::Vector4f min_p, max_p;
    // Get the minimum and maximum dimensions
    if (!filter_field_name_.empty())  // If we don't want to process the entire cloud...
        pcl::getMinMax3D<PointT>(input_, filter_field_name_, static_cast<float>(filter_limit_min_),
                                 static_cast<float>(filter_limit_max_), min_p, max_p,
                                 filter_limit_negative_);
    else
        pcl::getMinMax3D<PointT>(*input_, min_p, max_p);

    // Check that the leaf size is not too small, given the size of the data
    int64_t dx = static_cast<int64_t>((max_p[0] - min_p[0]) * inverse_leaf_size_[0]) + 1;
    int64_t dy = static_cast<int64_t>((max_p[1] - min_p[1]) * inverse_leaf_size_[1]) + 1;
    int64_t dz = static_cast<int64_t>((max_p[2] - min_p[2]) * inverse_leaf_size_[2]) + 1;

    if ((dx * dy * dz) > std::numeric_limits<int32_t>::max()) {
        PCL_WARN(
            "[pcl::%s::applyFilter] Leaf size is too small for the input dataset. Integer indices "
            "would overflow.",
            getClassName().c_str());
        output.clear();
        return;
    }

    // Compute the minimum and maximum bounding box values
    min_b_[0] = static_cast<int>(floor(min_p[0] * inverse_leaf_size_[0]));
    max_b_[0] = static_cast<int>(floor(max_p[0] * inverse_leaf_size_[0]));
    min_b_[1] = static_cast<int>(floor(min_p[1] * inverse_leaf_size_[1]));
    max_b_[1] = static_cast<int>(floor(max_p[1] * inverse_leaf_size_[1]));
    min_b_[2] = static_cast<int>(floor(min_p[2] * inverse_leaf_size_[2]));
    max_b_[2] = static_cast<int>(floor(max_p[2] * inverse_leaf_size_[2]));

    // Compute the number of divisions needed along all axis
    div_b_ = max_b_ - min_b_ + Eigen::Vector4i::Ones();
    div_b_[3] = 0;

    // Clear the leaves
    leaves_.clear();
    //  leaves_.reserve(8192);

    // Set up the division multiplier
    divb_mul_ = Eigen::Vector4i(1, div_b_[0], div_b_[0] * div_b_[1], 0);

    int centroid_size = 4;

    if (downsample_all_data_) centroid_size = boost::mpl::size<FieldList>::value;

    // ---[ RGB special case
    std::vector<pcl::PCLPointField> fields;
    int rgba_index = -1;
    rgba_index = pcl::getFieldIndex<PointT>("rgb", fields);
    if (rgba_index == -1) rgba_index = pcl::getFieldIndex<PointT>("rgba", fields);
    if (rgba_index >= 0) {
        rgba_index = fields[rgba_index].offset;
        centroid_size += 3;
    }

    // If we don't want to process the entire cloud, but rather filter points far away from the
    // viewpoint first...
    if (!filter_field_name_.empty()) {
        // Get the distance field index
        std::vector<pcl::PCLPointField> fields;
        int distance_idx = pcl::getFieldIndex<PointT>(filter_field_name_, fields);
        if (distance_idx == -1)
            PCL_WARN("[pcl::%s::applyFilter] Invalid filter field name. Index is %d.\n",
                     getClassName().c_str(), distance_idx);

        // First pass: go over all points and insert them into the right leaf
        for (size_t cp = 0; cp < input_->points.size(); ++cp) {
            if (!input_->is_dense)
                // Check if the point is invalid
                if (!std::isfinite(input_->points[cp].x) || !std::isfinite(input_->points[cp].y) ||
                    !std::isfinite(input_->points[cp].z))
                    continue;

            // Get the distance value
            const uint8_t *pt_data = reinterpret_cast<const uint8_t *>(&input_->points[cp]);
            float distance_value = 0;
            memcpy(&distance_value, pt_data + fields[distance_idx].offset, sizeof(float));

            if (filter_limit_negative_) {
                // Use a threshold for cutting out points which inside the interval
                if ((distance_value < filter_limit_max_) && (distance_value > filter_limit_min_))
                    continue;
            } else {
                // Use a threshold for cutting out points which are too close/far away
                if ((distance_value > filter_limit_max_) || (distance_value < filter_limit_min_))
                    continue;
            }

            int ijk0 = static_cast<int>(floor(input_->points[cp].x * inverse_leaf_size_[0]) -
                                        static_cast<float>(min_b_[0]));
            int ijk1 = static_cast<int>(floor(input_->points[cp].y * inverse_leaf_size_[1]) -
                                        static_cast<float>(min_b_[1]));
            int ijk2 = static_cast<int>(floor(input_->points[cp].z * inverse_leaf_size_[2]) -
                                        static_cast<float>(min_b_[2]));

            // Compute the centroid leaf index
            int idx = ijk0 * divb_mul_[0] + ijk1 * divb_mul_[1] + ijk2 * divb_mul_[2];

            Leaf &leaf = leaves_[idx];
            if (leaf.nr_points == 0) {
                leaf.centroid.resize(centroid_size);
                leaf.centroid.setZero();
            }

            Eigen::Vector3d pt3d(input_->points[cp].x, input_->points[cp].y, input_->points[cp].z);
            // Accumulate point sum for centroid calculation
            leaf.mean_ += pt3d;
            // Accumulate x*xT for single pass covariance calculation
            leaf.cov_ += pt3d * pt3d.transpose();

            // Do we need to process all the fields?
            if (!downsample_all_data_) {
                Eigen::Vector4f pt(input_->points[cp].x, input_->points[cp].y, input_->points[cp].z,
                                   0);
                leaf.centroid.template head<4>() += pt;
            } else {
                // Copy all the fields
                Eigen::VectorXf centroid = Eigen::VectorXf::Zero(centroid_size);
                // ---[ RGB special case
                if (rgba_index >= 0) {
                    // fill r/g/b data
                    int rgb;
                    memcpy(&rgb, reinterpret_cast<const char *>(&input_->points[cp]) + rgba_index,
                           sizeof(int));
                    centroid[centroid_size - 3] = static_cast<float>((rgb >> 16) & 0x0000ff);
                    centroid[centroid_size - 2] = static_cast<float>((rgb >> 8) & 0x0000ff);
                    centroid[centroid_size - 1] = static_cast<float>((rgb) & 0x0000ff);
                }
                pcl::for_each_type<FieldList>(
                    pcl::NdCopyPointEigenFunctor<PointT>(input_->points[cp], centroid));
                leaf.centroid += centroid;
            }
            ++leaf.nr_points;

            leaf.pointList_.push_back(input_->points[cp]);
        }
    }
    // No distance filtering, process all data
    else {
        // First pass: go over all points and insert them into the right leaf
        for (size_t cp = 0; cp < input_->points.size(); ++cp) {
            if (!input_->is_dense)
                // Check if the point is invalid
                if (!std::isfinite(input_->points[cp].x) || !std::isfinite(input_->points[cp].y) ||
                    !std::isfinite(input_->points[cp].z))
                    continue;

            int ijk0 = static_cast<int>(floor(input_->points[cp].x * inverse_leaf_size_[0]) -
                                        static_cast<float>(min_b_[0]));
            int ijk1 = static_cast<int>(floor(input_->points[cp].y * inverse_leaf_size_[1]) -
                                        static_cast<float>(min_b_[1]));
            int ijk2 = static_cast<int>(floor(input_->points[cp].z * inverse_leaf_size_[2]) -
                                        static_cast<float>(min_b_[2]));

            // Compute the centroid leaf index
            int idx = ijk0 * divb_mul_[0] + ijk1 * divb_mul_[1] + ijk2 * divb_mul_[2];

            // int idx = (((input_->points[cp].getArray4fmap () * inverse_leaf_size_).template
            // cast<int> ()).matrix () - min_b_).dot (divb_mul_);
            Leaf &leaf = leaves_[idx];
            if (leaf.nr_points == 0) {
                leaf.centroid.resize(centroid_size);
                leaf.centroid.setZero();
            }

            Eigen::Vector3d pt3d(input_->points[cp].x, input_->points[cp].y, input_->points[cp].z);
            // Accumulate point sum for centroid calculation
            leaf.mean_ += pt3d;
            // Accumulate x*xT for single pass covariance calculation
            leaf.cov_ += pt3d * pt3d.transpose();

            // Do we need to process all the fields?
            if (!downsample_all_data_) {
                Eigen::Vector4f pt(input_->points[cp].x, input_->points[cp].y, input_->points[cp].z,
                                   0);
                leaf.centroid.template head<4>() += pt;
            } else {
                // Copy all the fields
                Eigen::VectorXf centroid = Eigen::VectorXf::Zero(centroid_size);
                // ---[ RGB special case
                if (rgba_index >= 0) {
                    // Fill r/g/b data, assuming that the order is BGRA
                    int rgb;
                    memcpy(&rgb, reinterpret_cast<const char *>(&input_->points[cp]) + rgba_index,
                           sizeof(int));
                    centroid[centroid_size - 3] = static_cast<float>((rgb >> 16) & 0x0000ff);
                    centroid[centroid_size - 2] = static_cast<float>((rgb >> 8) & 0x0000ff);
                    centroid[centroid_size - 1] = static_cast<float>((rgb) & 0x0000ff);
                }
                pcl::for_each_type<FieldList>(
                    pcl::NdCopyPointEigenFunctor<PointT>(input_->points[cp], centroid));
                leaf.centroid += centroid;
            }
            ++leaf.nr_points;

            leaf.pointList_.push_back(input_->points[cp]);
        }
    }

    // Second pass: go over all leaves and compute centroids and covariance matrices
    output.points.reserve(leaves_.size());
    if (searchable_) voxel_centroids_leaf_indices_.reserve(leaves_.size());
    int cp = 0;
    if (save_leaf_layout_) leaf_layout_.resize(div_b_[0] * div_b_[1] * div_b_[2], -1);

    // Eigen values and vectors calculated to prevent near singluar matrices
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver;
    Eigen::Matrix3d eigen_val;
    Eigen::Vector3d pt_sum;

    // Eigen values less than a threshold of max eigen value are inflated to a set fraction of the
    // max eigen value.
    double min_covar_eigvalue;

    for (auto it = leaves_.begin(); it != leaves_.end(); ++it) {
        // Normalize the centroid
        Leaf &leaf = it->second;

        // Normalize the centroid
        leaf.centroid /= static_cast<float>(leaf.nr_points);
        // Point sum used for single pass covariance calculation
        pt_sum = leaf.mean_;
        // Normalize mean
        leaf.mean_ /= leaf.nr_points;

        // If the voxel contains sufficient points, its covariance is calculated and is added to the
        // voxel centroids and output clouds. Points with less than the minimum points will have a
        // can not be accuratly approximated using a normal distribution.
        if (leaf.nr_points >= min_points_per_voxel_) {
            if (save_leaf_layout_) leaf_layout_[it->first] = cp++;

            output.push_back(PointT());

            // Do we need to process all the fields?
            if (!downsample_all_data_) {
                output.points.back().x = leaf.centroid[0];
                output.points.back().y = leaf.centroid[1];
                output.points.back().z = leaf.centroid[2];
            } else {
                pcl::for_each_type<FieldList>(
                    pcl::NdCopyEigenPointFunctor<PointT>(leaf.centroid, output.back()));
                // ---[ RGB special case
                if (rgba_index >= 0) {
                    // pack r/g/b into rgb
                    float r = leaf.centroid[centroid_size - 3],
                          g = leaf.centroid[centroid_size - 2],
                          b = leaf.centroid[centroid_size - 1];
                    int rgb = (static_cast<int>(r)) << 16 | (static_cast<int>(g)) << 8 |
                              (static_cast<int>(b));
                    memcpy(reinterpret_cast<char *>(&output.points.back()) + rgba_index, &rgb,
                           sizeof(float));
                }
            }

            // Stores the voxel indice for fast access searching
            if (searchable_) voxel_centroids_leaf_indices_.push_back(static_cast<int>(it->first));

            // Single pass covariance calculation
            leaf.cov_ = (leaf.cov_ - 2 * (pt_sum * leaf.mean_.transpose())) / leaf.nr_points +
                        leaf.mean_ * leaf.mean_.transpose();
            leaf.cov_ *= (leaf.nr_points - 1.0) / leaf.nr_points;

            // Normalize Eigen Val such that max no more than 100x min.
            eigensolver.compute(leaf.cov_);
            eigen_val = eigensolver.eigenvalues().asDiagonal();
            leaf.evecs_ = eigensolver.eigenvectors();

            if (eigen_val(0, 0) < 0 || eigen_val(1, 1) < 0 || eigen_val(2, 2) <= 0) {
                leaf.nr_points = -1;
                continue;
            }

            // Avoids matrices near singularities (eq 6.11)[Magnusson 2009]

            min_covar_eigvalue = min_covar_eigvalue_mult_ * eigen_val(2, 2);
            if (eigen_val(0, 0) < min_covar_eigvalue) {
                eigen_val(0, 0) = min_covar_eigvalue;

                if (eigen_val(1, 1) < min_covar_eigvalue) {
                    eigen_val(1, 1) = min_covar_eigvalue;
                }

                leaf.cov_ = leaf.evecs_ * eigen_val * leaf.evecs_.inverse();
            }
            leaf.evals_ = eigen_val.diagonal();

            leaf.icov_ = leaf.cov_.inverse();
            if (leaf.icov_.maxCoeff() == std::numeric_limits<float>::infinity() ||
                leaf.icov_.minCoeff() == -std::numeric_limits<float>::infinity()) {
                leaf.nr_points = -1;
            }
        }
    }

    output.width = static_cast<uint32_t>(output.points.size());
}

//////////////////////////////////////////////////////////////////////////////////////////
template <typename PointT>
int pclomp::VoxelGridCovariance<PointT>::getNeighborhoodAtPoint(
    const Eigen::MatrixXi &relative_coordinates,
    const PointT &reference_point,
    std::vector<LeafConstPtr> &neighbors) const {
    neighbors.clear();

    // Find displacement coordinates
    Eigen::Vector4i ijk(static_cast<int>(floor(reference_point.x / leaf_size_[0])),
                        static_cast<int>(floor(reference_point.y / leaf_size_[1])),
                        static_cast<int>(floor(reference_point.z / leaf_size_[2])), 0);
    Eigen::Array4i diff2min = min_b_ - ijk;
    Eigen::Array4i diff2max = max_b_ - ijk;
    neighbors.reserve(relative_coordinates.cols());

    // Check each neighbor to see if it is occupied and contains sufficient points
    // Slower than radius search because needs to check 26 indices
    for (int ni = 0; ni < relative_coordinates.cols(); ni++) {
        Eigen::Vector4i displacement =
            (Eigen::Vector4i() << relative_coordinates.col(ni), 0).finished();
        // Checking if the specified cell is in the grid
        if ((diff2min <= displacement.array()).all() && (diff2max >= displacement.array()).all()) {
            auto leaf_iter = leaves_.find(((ijk + displacement - min_b_).dot(divb_mul_)));
            if (leaf_iter != leaves_.end() &&
                leaf_iter->second.nr_points >= min_points_per_voxel_) {
                LeafConstPtr leaf = &(leaf_iter->second);
                neighbors.push_back(leaf);
            }
        }
    }

    return (static_cast<int>(neighbors.size()));
}

//////////////////////////////////////////////////////////////////////////////////////////
template <typename PointT>
int pclomp::VoxelGridCovariance<PointT>::getNeighborhoodAtPoint(
    const PointT &reference_point, std::vector<LeafConstPtr> &neighbors) const {
    neighbors.clear();

    // Find displacement coordinates
    Eigen::MatrixXi relative_coordinates = pcl::getAllNeighborCellIndices();
    return getNeighborhoodAtPoint(relative_coordinates, reference_point, neighbors);
}

//////////////////////////////////////////////////////////////////////////////////////////
template <typename PointT>
int pclomp::VoxelGridCovariance<PointT>::getNeighborhoodAtPoint7(
    const PointT &reference_point, std::vector<LeafConstPtr> &neighbors) const {
    neighbors.clear();

    Eigen::MatrixXi relative_coordinates(3, 7);
    relative_coordinates.setZero();
    relative_coordinates(0, 1) = 1;
    relative_coordinates(0, 2) = -1;
    relative_coordinates(1, 3) = 1;
    relative_coordinates(1, 4) = -1;
    relative_coordinates(2, 5) = 1;
    relative_coordinates(2, 6) = -1;

    return getNeighborhoodAtPoint(relative_coordinates, reference_point, neighbors);
}

//////////////////////////////////////////////////////////////////////////////////////////
template <typename PointT>
int pclomp::VoxelGridCovariance<PointT>::getNeighborhoodAtPoint1(
    const PointT &reference_point, std::vector<LeafConstPtr> &neighbors) const {
    neighbors.clear();
    return getNeighborhoodAtPoint(Eigen::MatrixXi::Zero(3, 1), reference_point, neighbors);
}

//////////////////////////////////////////////////////////////////////////////////////////
template <typename PointT>
void pclomp::VoxelGridCovariance<PointT>::getDisplayCloud(
    pcl::PointCloud<pcl::PointXYZ> &cell_cloud) {
    cell_cloud.clear();

    int pnt_per_cell = 1000;
    boost::mt19937 rng;
    boost::normal_distribution<> nd(0.0, leaf_size_.head(3).norm());
    boost::variate_generator<boost::mt19937 &, boost::normal_distribution<>> var_nor(rng, nd);

    Eigen::LLT<Eigen::Matrix3d> llt_of_cov;
    Eigen::Matrix3d cholesky_decomp;
    Eigen::Vector3d cell_mean;
    Eigen::Vector3d rand_point;
    Eigen::Vector3d dist_point;

    // Generate points for each occupied voxel with sufficient points.
    for (auto it = leaves_.begin(); it != leaves_.end(); ++it) {
        Leaf &leaf = it->second;

        if (leaf.nr_points >= min_points_per_voxel_) {
            cell_mean = leaf.mean_;
            llt_of_cov.compute(leaf.cov_);
            cholesky_decomp = llt_of_cov.matrixL();

            // Random points generated by sampling the normal distribution given by voxel mean and
            // covariance matrix
            for (int i = 0; i < pnt_per_cell; i++) {
                rand_point = Eigen::Vector3d(var_nor(), var_nor(), var_nor());
                dist_point = cell_mean + cholesky_decomp * rand_point;
                cell_cloud.push_back(pcl::PointXYZ(static_cast<float>(dist_point(0)),
                                                   static_cast<float>(dist_point(1)),
                                                   static_cast<float>(dist_point(2))));
            }
        }
    }
}

}  // namespace pclomp

#endif  // PCL_VOXEL_GRID_COVARIANCE_OMP_H_
