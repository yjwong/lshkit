SET(TOOLS fitdata lsh-run mplsh-run mplsh-predict mplsh-tune sketch-run forest-run scan embed txt2bin dump-query)
FOREACH(TOOL ${TOOLS})
ADD_EXECUTABLE(${TOOL} ${TOOL}.cpp)
TARGET_LINK_LIBRARIES(${TOOL} lshkit ${Boost_LIBRARIES} ${GSL_LIBRARIES})
ENDFOREACH(TOOL)

