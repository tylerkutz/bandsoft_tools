find_package(MySQL REQUIRED)
add_definitions(-DRCDB_MYSQL)

set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/hipolib ${MYSQL_INCLUDE_DIR})

include_directories(${INCLUDES})

set(CMAKE_VERBOSE_MAKEFILE ON)
add_executable(neutrons neutrons.cpp readhipo_helper.cpp )
target_link_libraries(neutrons ${ROOT_LIBRARIES} Hipo4 Banks4 bandhit_tree genpart_tree -lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )
install( TARGETS neutrons
		DESTINATION ${CMAKE_SOURCE_DIR}/bin )

add_executable(electrons electrons.cpp readhipo_helper.cpp )
target_link_libraries(electrons ${ROOT_LIBRARIES} Hipo4 Banks4 clashit_tree genpart_tree clas12fiducial -lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )
install( TARGETS electrons
		DESTINATION ${CMAKE_SOURCE_DIR}/bin )

add_executable(gen_electrons gen_electrons.cpp)
target_link_libraries(gen_electrons ${ROOT_LIBRARIES} Hipo4 Banks4 clashit_tree genpart_tree clas12fiducial -lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )
install( TARGETS gen_electrons
		DESTINATION ${CMAKE_SOURCE_DIR}/bin )

add_executable(background background.cpp)
target_link_libraries(background ${ROOT_LIBRARIES} Hipo4 Banks4 clashit_tree genpart_tree clas12fiducial -lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )
install( TARGETS background
		DESTINATION ${CMAKE_SOURCE_DIR}/bin )


#message("Searching for Libraries")
#find_library(Test clashit_tree)
#find_library(Test clashit_tree HINTS /w/hallb-scifs17exp/clas12/users/awild/BAND_Project/Build/clashit)
#message(${Test})
#if(NOT Test)
#  message(FATAL_ERROR)
#endif()


add_executable(tagged_dis tagged_dis.cpp readhipo_helper.cpp )
target_link_libraries(tagged_dis ${ROOT_LIBRARIES} Hipo4 Banks4  bandhit_tree clashit_tree taghit_tree genpart_tree -lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )
install( TARGETS tagged_dis
	DESTINATION ${CMAKE_SOURCE_DIR}/bin )

add_executable(ncalibration ncalibration.cpp readhipo_helper.cpp )
target_link_libraries(ncalibration ${ROOT_LIBRARIES} Hipo4 Banks4  bandhit_tree clashit_tree  -lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )
install( TARGETS ncalibration
	DESTINATION ${CMAKE_SOURCE_DIR}/bin )

#add_executable(EC_cuts EC_cuts.cpp readhipo_helper.cpp )
#target_link_libraries(EC_cuts ${ROOT_LIBRARIES} Hipo4 Banks4 clashit_tree clas12fiducial e_pid -lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )
#install( TARGETS EC_cuts
#		DESTINATION ${CMAKE_SOURCE_DIR}/bin )

#add_executable(SF_cuts_fit SF_cuts_fit.cpp readhipo_helper.cpp )
#target_link_libraries(SF_cuts_fit ${ROOT_LIBRARIES} Hipo4 Banks4 clashit_tree clas12fiducial -lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )
#install( TARGETS SF_cuts_fit
#		DESTINATION ${CMAKE_SOURCE_DIR}/bin )

#add_executable(SF_cuts_fit SF_cuts_fit.cpp readhipo_helper.cpp )
#target_link_libraries(SF_cuts_fit ${ROOT_LIBRARIES} Hipo4 Banks4 clashit_tree clas12fiducial -lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )
#install( TARGETS SF_cuts_fit
#		DESTINATION ${CMAKE_SOURCE_DIR}/bin )


add_executable(quasielastic quasielastic.cpp readhipo_helper.cpp )
target_link_libraries(quasielastic ${ROOT_LIBRARIES} Hipo4 Banks4  bandhit_tree clashit_tree -lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )
install( TARGETS quasielastic
		DESTINATION ${CMAKE_SOURCE_DIR}/bin )
		
add_executable(DC_fid_cuts DC_fid_cuts.cpp)
target_link_libraries(DC_fid_cuts ${ROOT_LIBRARIES} Hipo4 Banks4 clashit_tree clas12fiducial DC_fiducial -lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )
install( TARGETS DC_fid_cuts
		DESTINATION ${CMAKE_SOURCE_DIR}/bin )
