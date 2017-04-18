#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>



//#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/Polyhedron_items_3.h>

#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>

//Sizing field
#include <CGAL/Mesh_3/experimental/Lipschitz_sizing_polyhedron.h>

// IO
#include <CGAL/IO/Polyhedron_iostream.h>
#include <list>
#include <string>
#include <CGAL/Mesh_3/Dump_c3t3.h>
#include <CGAL/IO/File_medit.h>
// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Mesh_polyhedron_3<K>::type Polyhedron;
typedef CGAL::Polyhedral_mesh_domain_with_features_3<K> Mesh_domain;
typedef K::FT FT;
#ifdef CGAL_CONCURRENT_MESH_3
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain,CGAL::Default,Concurrency_tag>::type Tr;


typedef CGAL::Mesh_complex_3_in_triangulation_3<
  Tr,Mesh_domain::Corner_index,Mesh_domain::Curve_segment_index> C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;




typedef CGAL::Mesh_3::Lipschitz_sizing<K, Mesh_domain, Mesh_domain::AABB_tree> Lip_sizing;


using namespace CGAL::parameters;


int main(int argc, char*argv[])
{
  const char* fname;
  const char* OutPath;
  const char* fname_b=argv[1];
  // Create input polyhedron from asc 
  std::list<Polyhedron*> poly_pointer_list; 
   
  
  Polyhedron bounding_polyhedron;
  std::ifstream input_b(fname_b);
  input_b >> bounding_polyhedron; 
  if(input_b.fail())
  {
     std::cerr << "Error: Cannot read file 1" <<  fname_b << std::endl;
     return EXIT_FAILURE;
  }
  input_b.close();  


  for (int i = 2; i < argc-1; i++) 
  { 
     fname = argv[i];
     std::ifstream input(fname);
     Polyhedron *polyhedron = new Polyhedron();
     input >> *polyhedron;
     if(input.fail())
     {
       std::cerr << "Error: Cannot read file 2" <<  fname << std::endl;
       return EXIT_FAILURE;
     }
     input.close();
     std::cout << polyhedron << std::endl;
     poly_pointer_list.push_back(polyhedron);

    
     
  }

  

  OutPath = argv[argc-1];

 
  

  // Create domain


  Mesh_domain domain(poly_pointer_list.begin(),poly_pointer_list.end(),bounding_polyhedron);	

  domain.make_surface_index();
  domain.detect_features(FT(30)) ;
  //domain.detect_borders();
  // not sure what this does
  Lip_sizing lip_sizing(domain, &domain.aabb_tree());
  FT min_size = 0.02;
  lip_sizing.add_parameters_for_subdomain(1,       //subdomain id
                                          0.3,     //k
                                          min_size,//min_size
                                          0.5);    //max_size

  // Mesh criteria (no cell_size set)
  Mesh_criteria criteria(facet_angle=25, facet_size=0.5, facet_distance=0.05,cell_radius_edge_ratio=3);

  
  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,perturb(),exude()); // perturb and exude ?? 

  // Output
 
  dump_c3t3(c3t3, "dump.mesh");
  std::ofstream medit_file(OutPath);
  //c3t3.output_to_medit(medit_file);
  output_to_medit(medit_file,c3t3,true,true);
  medit_file.close();

  return 0;
}
