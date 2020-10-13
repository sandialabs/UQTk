#ifndef PATHTOOLS_H
#define PATHTOOLS_H

#include <string>

namespace muq
{

namespace Utilities
{
    /**
       @brief Get the parent folder of a path
       @details Given a path of the form /a/b/c, return the parent path /a/b
\code
std::string base   = "/some/path";
std::string parent = GetParent(base);
\endcode
Parent now holds "/some"
     */	
    std::string GetParentPath(std::string const& base);


    /** 
	@brief Splits a string on the first forward slash.
	@details Given a path, like "/a/b/c", this function returns two strings, "/a", and "/b/c".
	@params[in] path The original plan to split.
	@return A pair of std::string holind the two components of the path.
    */
    std::pair<std::string, std::string> SplitString(std::string const& path);

} // namespace Utilities
} // namespace muq





#endif // #ifndef PATHTOOLS_H
