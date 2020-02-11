/*=========================================================================

  idpRegistrationUtilities.cxx

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/
#ifndef __idpRegistrationUtilities_cxx
#define __idpRegistrationUtilities_cxx

#include "idpRegistrationUtilities.h"


namespace itk
{
namespace idp
{

/*--------------------------
 * Helper function to change path between cluster and windows machines
 * --------------------------
 */
void ChangePaths( std::string & str )
{
#ifdef WIN32
  if (!itksys::SystemTools::StringStartsWith( str.c_str(), "//titan/CNS/" ) )
    {
    if ( str.find("/projects/devmouse/vol1/") != std::string::npos )
      {
      itksys::SystemTools::ReplaceString( str, "/projects/devmouse/vol1/", "//titan/CNS/devmouse/" );
      }
    else if ( str.find("/projects/0378/vol1/") != std::string::npos )
      {
      itksys::SystemTools::ReplaceString( str, "/projects/0378/vol1/", "//titan/CNS/0378/" ); 
      }
    else if ( str.find("/projects/0310/vol1/") != std::string::npos )
      {
      itksys::SystemTools::ReplaceString( str, "/projects/0310/vol1/", "//titan/CNS/0310/" );
      }
    else if ( str.find("/projects/aibssan/production30/") != std::string::npos )
      {
      itksys::SystemTools::ReplaceString( str, "/projects/aibssan/production30/", "//AIBS-WST-ISI-lb/ifs/data/aibssan/production30/" );  
      }
    else if ( str.find("/projects/aibssan/") != std::string::npos )
      {
      itksys::SystemTools::ReplaceString( str, "/projects/aibssan/", "//titan/CNS/" );  
      }
    else if ( str.find("/aibssan/") != std::string::npos )
      {
      itksys::SystemTools::ReplaceString( str, "/aibssan/", "//AIBS-WST-ISI-lb/ifs/data/aibssan/" );
      }
    }
#endif
}



Relation::Relation(double w,int v1,int v2){
	weight=w;
	vert1=v1;
	vert2=v2;
}

bool Relation::operator<(const Relation & b) const{
	return (this->weight>b.weight);
}

int newRelation(Relation edge1,Relation& edgeout,bool* vertices){
	bool new1=vertices[edge1.vert1]; 
	bool new2=vertices[edge1.vert2];
	int out1;
	if(new1^new2){
		if(new1){
			out1=edge1.vert2;
			edgeout=Relation(edge1.weight,edge1.vert1,edge1.vert2);
		}
		else {
			out1=edge1.vert1;
			edgeout=Relation(edge1.weight,edge1.vert2,edge1.vert1);
		}
	}
	else{
		out1=-1;
	}
	return out1;
}

} // end namespace idp
} //end namespace itk

#endif

  
