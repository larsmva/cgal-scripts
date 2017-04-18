#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>


#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Mesh_polyhedron_3.h>
#include <CGAL/Polyhedral_mesh_domain_3.h>
//#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>

//Sizing field

// IO
#include <CGAL/IO/Polyhedron_iostream.h>
#include <list>
#include <string>
#include <CGAL/Mesh_3/Dump_c3t3.h>
#include <CGAL/IO/File_medit.h>


// Domain Mesh_polyhedron_3
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, K> Polyhedral_mesh_domain_3;

// THIS IS A COPY OF A CGAL SCRIPT 
// wrapping polyhedral_mesh_domain_3, which has the class function is_in_domain, which is called on  
// within Labeled_mesh_domain_3. IS similar construction is possible for Polyhedral_mesh_domain_with_features 3?
// Possibility of a subclass of Polyhedral_mesh_domain_with_features 3 with different is_in_domain function.
namespace CGAL { 
        template<class Function_, class BGT> 
        class Polyhedral_vector_to_labeled_function_wrapper 
        { 
            public: 
                // Types 
                typedef int	 return_type; 
                typedef std::vector<Function_*>   Function_vector; 
                typedef typename BGT::Point_3	 Point_3; 
                /// Constructor 
                
                Polyhedral_vector_to_labeled_function_wrapper(std::vector<Function_*>& v) : function_vector_(v) {} 
                /// Destructor 
                ~Polyhedral_vector_to_labeled_function_wrapper() {} 

                /// Operator () 
                return_type operator()(const Point_3& p, const bool = true) const 
                { 
                    int nb_func = function_vector_.size(); 
                    if ( nb_func > 8 ) 
                    { 
                        CGAL_error_msg("We support at most 8 functions !"); 
                    } 

                    char bits = 0; 
                    for ( int i = 0 ; i < nb_func ; ++i ) 
                    {    
                        bits |= (( function_vector_[i]->is_in_domain_object()(p) > 0) << i ); 
                    } 
                    return ( static_cast<return_type>(bits) ); 
                } 
            private: 
                /// Functions to wrap 
                Function_vector function_vector_; 
        }; 
} 


#ifdef CGAL_CONCURRENT_MESH_3
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif
typedef CGAL::Polyhedral_vector_to_labeled_function_wrapper<Polyhedral_mesh_domain_3, K> Function_wrapper;
typedef Function_wrapper::Function_vector Function_vector; 
typedef CGAL::Labeled_mesh_domain_3<Function_wrapper, K> Mesh_domain; 

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain,K,Concurrency_tag>::type Tr;

typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

using namespace CGAL::parameters;


int main(int argc, char*argv[])
{
  const char* fname;
  const char* OutPath;
  const char* fname_b=argv[1];

  Function_vector v ; 
 
  Polyhedron bounding_polyhedron;
  std::ifstream input_b(fname_b);
  input_b >> bounding_polyhedron; 
 
  Polyhedral_mesh_domain_3 bounding_polyhedral_domain( bounding_polyhedron);




  if(input_b.fail())
  {
     std::cerr << "Error: Cannot read file 1" <<  fname_b << std::endl;
     return EXIT_FAILURE;
  }
  input_b.close();  

  v.push_back(&bounding_polyhedral_domain); 

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

     Polyhedral_mesh_domain_3 *polyhedral_domain = new Polyhedral_mesh_domain_3(*polyhedron);

     v.push_back(polyhedral_domain);
  
     //poly_pointer_list.push_back(polyhedron);

    
     
  }

  

  OutPath = argv[argc-1];


  // Create domain


  //Polyhedral_mesh_domain_3 p_domain(poly_pointer_list.begin(),poly_pointer_list.end(),bounding_polyhedron);	



  Mesh_domain domain(Function_wrapper(v),bounding_polyhedral_domain.bbox() );

  Mesh_criteria criteria(facet_angle=25, facet_size=0.5, facet_distance=0.1,cell_radius_edge_ratio=3);

  
  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria); // perturb and exude ?? 
  
  // Output
  dump_c3t3(c3t3, "red");
  //Maybe delete poly_point 
  std::ofstream *medit_file= new std::ofstream(OutPath);
  //c3t3.output_to_medit(medit_file);
  output_to_medit(*medit_file,c3t3,true,true);
  medit_file->close();

  return 0;
}
