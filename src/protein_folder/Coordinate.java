/*______________________________________________________________________________
|                                                                              |
|o------------------------------------>X<-------------------------------------o|
|                                                                              |
| ___________________________________________________________________________  |
| |   This program was written by:           4/28/2014                       | |
| |    Samuel David Lotz                                                     | |
| |     Title:                                                               | |
| |	 Protein_Folder                                                          | |
| |       Protein Folder is a project written under the guidance of
| |   Dr. Sam Thangiah at Slippery Rock University for a Machine learning 
| |   independent instruction course.  Code design and pseudo-
| |   code taken from Alena Shmygelska's dissertation "Novel Heuristic
| |   Search Methods for Protein Folding and Identification of Folding 
| |   Pathways" Sept 2006 University of British Columbia. The goal of this
| |   project was to implement machine learning strategies to the 
| |   biological problem of protein folding and to compare them to non-
| |   Learning stategies.
| | 
| |__________________________________________________________________________| |
|______________________________________________________________________________|*/
/* -------  Coordinate Class-------------------------------------------------------

	  The coordinate class is an encapsulated coordinate system.  It is encapsulated
	 so that it can be modified later to allow for multiple coordinate systems such 
	 as 3D.  Polymorphism and inheritance may come in handy, this is not an abstract 
	 class however.
	
	ATTRIBUTES:
		Coordinate variables:
			x : int -- the x-axis coordinate
			y : int -- the y-axis coordinate
			type : ResidueType -- the type of residue at this coordinate
				useful for printing
	
	CONSTRUCTORS:
	Coordinate() : No-argument -- creates a coordinate with x=y=0 and an uninitialized
		residue type
	Coodinate(int x, int y) : coordinate specified constructor -- input parameter sets the 
		x and y coordinates and does not initialize the residue type.
	
	
	METHODS:
	GETTERS
	SETTERS
	equals(Coordinate comparedCoordinate): boolean -- compares the calling coordinate to
		the one specified in the parameters.  If they are at the same location will return true
		regardless of type. False if they are not at the sae location.
	incrementers 
		incY, incX -- increase the coordinate by one LATTICE_UNIT defined in Param
		decY, decX -- decrease the coodinate by one LATTICE_UNIT defined in Param
	stringCoord(): String -- returns a formatted string of this coordinate
		e.g.   x,y (type)
	
	
	
	Version 1.0 beta 
		- Only 2D proteins are supported.  Coordinates are x and y
  ------------------------------------------------------------------------------
 */

package protein_folder;

import protein_folder.Param.ResidueType;

public class Coordinate{
	
//==================ATTRIBUTES======================================================
	private int x;
	private int y;
	//and the residue type to make calculations easy
	private ResidueType type;
//==================================================================================
	
//////__________Constructors_________________________________________________________
	
//----No-Arg Constructor-----------------------------------------------------------
	
	public Coordinate(){
		x=0;
		y=0;
	}
//----End No-Arg Constructor-------------------------------------------------------

//----Coordinate Specified Constructor----------------------------------------------

	public Coordinate(int inX, int inY){
		x=inX;
		y=inY;
	}
//----End Coordinate Specified Constructor-------------------------------------------------------
	
//----Copy Constructor-----------------------------------------------------------
	//Doesnt work for some reason....
	/*copy consructor all coordinate values from another coordinate object
		public Coordinate(Coordinate inCoord){
			Coordinate outCoord = new Coordinate();
			outCoord.setX(inCoord.getX());
			outCoord.setY(inCoord.getY());
	}
	*/ 
//----End Copy Constructor-------------------------------------------------------

//////++++++++GETTERS & SETTERS+++++++++++++++++++++++++++++++++++++++++++++++++++++
	public int getX(){
		return x;
	}
	public int getY(){
		return y;
	}
	public ResidueType getType(){
		return type;
	}
	public void setX(int inX){
		x = inX;
	}
	public void setY(int inY){
		y = inY;
	}
	public void setType(ResidueType inType){
		this.type = inType;
	}
	
//////++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	
//==================METHODS======================================================
	
//----increments-----------------------------------------------------------
	public void incX(){
		x = x + Param.LATTICE_UNIT;
	}
	public void decX(){
		x = x - Param.LATTICE_UNIT;
	}
	public void incY(){
		y = y + Param.LATTICE_UNIT;
	}
	public void decY(){
		y = y - Param.LATTICE_UNIT;
	}
//----End Increments-------------------------------------------------------

//----equals-----------------------------------------------------------
	//checks to see if this coordinate is the same as another
	public boolean equals(Coordinate rhsCoord){
		if(this.x == rhsCoord.getX() && this.y == rhsCoord.getY())
			return true;
		else 
			return false;
	}
//----End equals-------------------------------------------------------
	
//----stringCoord-----------------------------------------------------------
	//returns a formatted coordinate "x,y"
	public String stringCoord(){
		String typeString = "invalid";
		switch(this.type){ case H: typeString = "H";
						break;
					  case P: typeString = "P";
					    break;
					  default: typeString = "invalid Type";
					    break;
		}//end switch
		
		return x + "," + y + " (" + typeString + ")";
	}
//----End stringCoord-------------------------------------------------------
	
//==================================================================================
	
}//end coordinate class