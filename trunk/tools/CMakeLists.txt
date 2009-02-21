SET(TOOLS fitdata lsh-run mplsh-run mplsh-predict mplsh-tune sketch-run scan embed)
FOREACH(TOOL ${TOOLS})
ADD_EXECUTABLE(${TOOL} ${TOOL}.cpp)
TARGET_LINK_LIBRARIES(${TOOL} lshkit ${Boost_LIBRARIES} ${GSL_LIBRARIES})
ENDFOREACH(TOOL)

