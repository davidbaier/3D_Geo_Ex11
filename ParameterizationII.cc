#include "ParameterizationII.hh"
#include <ACG/Utils/StopWatch.hh>

void ParameterizationII::map_suface_boundary_to_circle(const Point& _origin, const double _radius) {
    mesh_.request_vertex_texcoords2D();
    for(auto vh : mesh_.vertices())
        mesh_.set_texcoord2D(vh, Vec2d(0., 0.));

    //copy mesh_ to tex_mesh_
    for(auto vh : mesh_.vertices())
        tex_mesh_.add_vertex(mesh_.point(vh));
    for(auto fh : mesh_.faces()) {
        std::vector<OpenMesh::VertexHandle> fh_vhs;
        for(auto fv_it = mesh_.fv_iter(fh); fv_it.is_valid(); ++fv_it)
            fh_vhs.push_back(*fv_it);

        tex_mesh_.add_face(fh_vhs);
    }

    //Finding the boundary
    OpenMesh::VertexHandle vh_start;
    for(auto vh : mesh_.vertices()) {
        if(mesh_.is_boundary(vh)) {
            //start from a manifold vertex and valence larger than 2
            if(mesh_.is_manifold(vh) && mesh_.valence(vh) > 2){
                vh_start = vh;
                break;
            }
        }
    }
    if(!vh_start.is_valid()) {
        std::cerr<<"\nProbably the input mesh has no boundary!";
        return;
    }

    //Start from a boundary vertex. Traverse all the boundary edges
    //and map the boundary vertices onto a circle with radius _radius,
    //the interior vertices to the circle center (0,0) in XY plane.
    // ------------- IMPLEMENT HERE ---------
    
    for(auto vh : mesh_.vertices())
        mesh_.set_texcoord2D(vh, Vec2d(_origin[0], _origin[1]));

    // all the commented out functions are there for an approach with unordered_map.
    
    //std::unordered_map<double, OpenMesh::VertexHandle> boundary_vertices;
    std::vector<OpenMesh::VertexHandle> boundary_vertices;
    std::vector<double> boundary_length;

    double total_length = 0.0;
    //boundary_vertices.insert(std::make_pair(0, vh_start));
    //boundary_vertices.emplace_back(vh_start);
    while(vh_start.is_valid()){
        for(auto v_hitr = mesh_.voh_iter(vh_start); mesh_.is_boundary(*v_hitr); v_hitr++){
            auto vert_handle = mesh_.to_vertex_handle(*v_hitr);
            double length = (mesh_.point(vert_handle) - mesh_.point(vh_start)).norm();
            if(std::find(boundary_vertices.begin(), boundary_vertices.end(), vert_handle) == boundary_vertices.end()){
            //if(boundary_vertices.find(length) == boundary_vertices.end()){
                //boundary_vertices.insert(std::make_pair(length, vert_handle));
                boundary_vertices.emplace_back(vert_handle);
                boundary_length.emplace_back(length);
                total_length += length;
                vh_start = vert_handle;
            }else{
                vh_start = OpenMesh::VertexHandle();
            }
        }
    }

    std::cout << _origin << std::endl;
    double angle_summed = 0.0;
    //for(auto boundary : boundary_vertices){
    for(int i = 0; i < boundary_vertices.size(); ++i){
        //auto length = boundary.first;
        //auto vertices = boundary.second;
        auto length = boundary_length[i];
        auto vertices = boundary_vertices[i];

        auto ratio = length/total_length*2.0*M_PI;
        auto angle = angle_summed + ratio;
        angle_summed += ratio;
        Vec2d coords (_origin[0]+_radius*std::cos(angle),_origin[1]+_radius*std::sin(angle));
        mesh_.set_texcoord2D(vertices, coords);
    }

    // ------------- IMPLEMENT HERE ---------

    //Update the texture mesh
    for (auto vh : mesh_.vertices()) {
        auto tex_coord = mesh_.texcoord2D(vh);
        tex_mesh_.set_point(vh, Point(tex_coord[0], tex_coord[1], 0) + _origin);
    }
}


// ======================================================================
// EXERCISE 1.2 Interactively smoothing the texture
// ========================================================================
void ParameterizationII::explicitly_smooth_texture(const Point& _origin, const int _num_iters)
{
    //See the handout
    // ------------- IMPLEMENT HERE ---------
    for (int iter=0; iter<_num_iters; ++iter) {
    calc_edges_weights();
        for(auto vh : mesh_.vertices()){
            if(!mesh_.is_boundary(vh)){
                auto total_weight = 0.0;
                Point summed_weigthed_neighbors = Point(0.);

                for(auto v_hitr = mesh_.voh_iter(vh); v_hitr.is_valid(); ++v_hitr){
                    auto edge_handle = mesh_.edge_handle(*v_hitr);
                    auto edge_weight = mesh_.property(edge_weight_, edge_handle);

                    auto neighbor_vh = mesh_.to_vertex_handle(*v_hitr);
                    if(edge_weight < 0.0)
                        std::cerr << "For vertex " << neighbor_vh << " the edge weight is smaller than 0.0" << std::endl;
                    total_weight += edge_weight; 

                    summed_weigthed_neighbors += edge_weight * mesh_.point(neighbor_vh);
                }

                if(total_weight > 1.0)
                    std::cerr << "For vertex " << vh << "the total weight is above 1.0" << std::endl;
                summed_weigthed_neighbors /= total_weight;
                mesh_.point(vh) = summed_weigthed_neighbors;
            }

        }
    }


    // ------------- IMPLEMENT HERE ---------

    //Update the texture mesh
    for (auto vh : mesh_.vertices()) {
        auto tex_coord = mesh_.texcoord2D(vh);
        tex_mesh_.set_point(vh, Point(tex_coord[0], tex_coord[1], 0) + _origin);
    }
}


// ======================================================================
// EXERCISE 1.3 Implicitly smoothing the texture
// ========================================================================
void ParameterizationII::implicitly_smooth_texture(const Point& _origin)
{
    //See the handout
    // ------------- IMPLEMENT HERE ---------


    // ------------- IMPLEMENT HERE ---------

    //Update the texture mesh
    for (auto vh : mesh_.vertices()) {
        auto tex_coord = mesh_.texcoord2D(vh);
        tex_mesh_.set_point(vh, Point(tex_coord[0], tex_coord[1], 0) + _origin);
    }
}


// ======================================================================
// EXERCISE 2 Minimal Surfaces
// ======================================================================
void ParameterizationII::minimal_surface() {
    //See the handout
    // ------------- IMPLEMENT HERE ---------


    // ------------- IMPLEMENT HERE ---------

    mesh_.update_normals();
}