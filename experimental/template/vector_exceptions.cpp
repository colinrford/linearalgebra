
module;

import <exception>;

export module la.vectors.exceptions;

export namespace la
{
namespace vectors
{
  struct vector_exception : public std::exception {

    struct non_pos : public std::exception {
      const char* what () const throw () {
        return "The dimension of the vector must be at least 1.\n";
      }
    };
    struct div_by_zero : public std::exception {
      const char* what () const throw () {
        return "Division by zero, ' 0 ', is undefined!\n";
      }
    };
    struct cross_undef : public std::exception {
      const char* what () const throw () {
        return "This cross product is not defined in dimensions other than 3. \
                Wedges later.\n";
      }
    };
    struct out_of_bounds : public std::exception {
      const char* what () const throw () {
        return "Index exceeds dimension of vector.\n";
      }
    };
    struct null_construction : public std::exception {
      const char* what () const throw () {
        return "Construction of a vector with nullptr is not allowed.\n";
      }
    };
    struct plus_equals_unequal_dim : public std::exception {
      const char* what () const throw () {
        return "operator+= is more restricted than others, it requires the \
                vectors have equal dimension.\n";
      }
    };
  };

} // end namespace vector
} // end namespace la
