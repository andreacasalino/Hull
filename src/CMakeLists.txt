set(PROJECT_SHORTNAME Hull)

MakeLibrary(${PROJECT_SHORTNAME} header)

target_link_libraries(${PROJECT_SHORTNAME} PUBLIC
${EigenName_To_Link}
)