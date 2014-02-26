.cc.o:
	$(CXX) $(ALL_CPPFLAGS) $(ALL_CXXFLAGS) -c $< -o $@
.cpp.o:
	$(CXX) $(ALL_CPPFLAGS) $(ALL_CXXFLAGS) -c $< -o $@
.cc.d:
	@echo "making $@ from $<"
	@$(CXX) -M $(ALL_CPPFLAGS) $< > $@
.cpp.d:
	@echo "making $@ from $<"
	@$(CXX) -M $(ALL_CPPFLAGS) $< > $@
