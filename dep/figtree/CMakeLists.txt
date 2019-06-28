

FILE(GLOB fgsrcs "*.cpp")
add_library(depfigtree ${fgsrcs})

include_directories (../ann)

INSTALL(TARGETS depfigtree DESTINATION lib)

# Install the header files
SET(figtree_HEADERS 
  figtree_internal.h
  figtree.h
  KCenterClustering.h
  )

INSTALL(FILES ${figtree_HEADERS} DESTINATION include/dep)

