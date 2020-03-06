

FILE(GLOB ansrcs "*.cpp")
add_library(depann ${ansrcs})


INSTALL(TARGETS depann DESTINATION lib)

# Install the header files
SET(ann_HEADERS 
  ANN.h
  ANNx.h
  ANNperf.h
  bd_tree.h
  kd_fix_rad_search.h
  kd_pr_search.h
  kd_search.h
  kd_split.h
  kd_tree.h
  kd_util.h
  pr_queue_k.h
  pr_queue.h
  )

INSTALL(FILES ${ann_HEADERS} DESTINATION include/dep)

