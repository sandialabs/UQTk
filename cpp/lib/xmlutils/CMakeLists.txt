
FILE(GLOB slsrc "*.cpp")

SET(xml_HEADERS
  XMLUtils.h
  XMLExpatParser.h
  XMLParser.h
  XMLElement.h
  Object.h
  RefPtr.h
  MyException.h
  XMLAttributeList.h
  )

add_library(uqtkxmlutils ${slsrc})

# Install the library
INSTALL(TARGETS uqtkxmlutils DESTINATION lib)

# Install the header files
INSTALL(FILES ${xml_HEADERS} DESTINATION include/uqtk)


