/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, development version     *
*                (c) 2006-2016 INRIA, USTL, UJF, CNRS, MGH                    *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU General Public License as published by the Free  *
* Software Foundation; either version 2 of the License, or (at your option)   *
* any later version.                                                          *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for    *
* more details.                                                               *
*                                                                             *
* You should have received a copy of the GNU General Public License along     *
* with this program; if not, write to the Free Software Foundation, Inc., 51  *
* Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.                   *
*******************************************************************************
*                            SOFA :: Applications                             *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/helper/OptionsGroup.h>
#include <sofa/simulation/AnimateBeginEvent.h>
#include <SofaTest/Sofa_test.h>
#include <sofa/core/ObjectFactory.h>


#include <SofaTest/TestMessageHandler.h>
using sofa::helper::logging::ExpectMessage ;
using sofa::helper::logging::Message ;

#include <SofaSimulationGraph/DAGSimulation.h>
using sofa::core::objectmodel::ComponentState ;
using sofa::core::objectmodel::BaseObject ;
using sofa::simulation::Simulation ;
using sofa::simulation::Node ;

#include <SofaSimulationCommon/SceneLoaderXML.h>
using sofa::simulation::SceneLoaderXML ;
using sofa::core::ExecParams ;

#include <sofa/helper/system/config.h>

#include <gtest/gtest.h>

#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>


//#ifdef SOFA_BUILD_FAKEPLUGIN
//#define SOFA_FAKEPlugin_API SOFA_EXPORT_DYNAMIC_LIBRARY
//#else
//#define SOFA_FAKEPlugin_API SOFA_IMPORT_DYNAMIC_LIBRARY
//#endif

//extern "C" {

//SOFA_FAKEPlugin_API void initExternalModule();

//SOFA_FAKEPlugin_API const char* getModuleName();

//SOFA_FAKEPlugin_API const char* getModuleVersion();

//SOFA_FAKEPlugin_API const char* getModuleLicense();

//SOFA_FAKEPlugin_API const char* getModuleDescription();

//SOFA_FAKEPlugin_API const char* getModuleComponentList();

//}

//extern "C" {

//void initExternalModule()
//{
//}

//const char* getModuleName()
//{
//		return "Fake plugin";
//}

//const char* getModuleVersion()
//{
//		return "1.0";
//}

//const char* getModuleLicense()
//{
//		return "None";
//}

//const char* getModuleDescription()
//{
//		return "";
//}

//const char* getModuleComponentList()
//{
//		return "TestObject";
//}

//}



namespace sofa {

class TestObject : public core::objectmodel::BaseObject
{

public:
	SOFA_CLASS(TestObject, core::objectmodel::BaseObject);

	Data<helper::OptionsGroup> opts;
	Data<std::string> optsValues;

	TestObject()
	: Inherit1(),
	opts(initData(&opts, "opt", "OptionsGroup initial value")),
	optsValues(initData(&optsValues, std::string("A;B;C;D;E;F"), "values", "All the potential values for the combobox (could imagine a directory path containg files, files would be the different values"))
	{
		std::cout << "lol" << std::endl;
	}

	~TestObject() {}

	template<typename Out>
	void split(const std::string &s, char delim, Out result) {
			std::stringstream ss;
			ss.str(s);
			std::string item;
			while (std::getline(ss, item, delim)) {
					*(result++) = item;
			}
	}
	std::vector<std::string> split(const std::string &s, char delim) {
			std::vector<std::string> elems;
			split(s, delim, std::back_inserter(elems));
			return elems;
	}

	void init()
	{
		// Sets the default value of the OptionsGroup, as defined in the scene file
		opts.beginWriteOnly()->setSelectedItemToDefault();
		opts.endEdit();

		// Assigns the values retrieved from the other Data<>
		std::vector<std::string> tokens = split(optsValues.getValue(), ';');

		helper::OptionsGroup* t = opts.beginEdit();
		t->setNbItems(unsigned(tokens.size()));
		unsigned i = 0;
		for (std::string& s : tokens)
			t->setItemName(i++, s);
		opts.endEdit();
	}
};

SOFA_DECL_CLASS(TestObject)

int TestObjectClass =
		core::RegisterObject("OptionsGroup TestObject")
				.add<TestObject>();


struct OptionsGroup_test: public BaseSofa_test
{
		void testOptionsGroup()
		{
			this->clearSceneGraph();

			// This is a RAII message.
			ExpectMessage error(Message::Error) ;

			std::stringstream scene ;
			scene << "<Node 	name='Root'>     \n"
							 "  <TestObject name='to' opt='C'/>  \n"
							 "</Node>                  \n" ;

			Node::SPtr root = SceneLoaderXML::loadFromMemory ("testscene",
																												scene.str().c_str(),
																												scene.str().size()) ;

			root->init(ExecParams::defaultInstance()) ;


			TestObject* testobject = dynamic_cast<TestObject*>(root->getObject("to"));
			EXPECT_NE(testobject, nullptr) ;

			/// The value C, the 3rd item in the 0-indexed OptionsGroup, is the default value in the scene and should be properly recovered,
			EXPECT_EQ(testobject->opts.getValue().getSelectedId(), 2) ;
		}
};

// Test
TEST_F(OptionsGroup_test, testOptionsGroup)
{
		this->testOptionsGroup();
}

} // namespace sofa
