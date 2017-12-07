
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
/* -------  Residue Class-------------------------------------------------------

	 This class is intended to encapsulate information associated with each 
	 residue in the sequence of a protein.  Such that a protein is a data
	 structure containing references to many residue objects.
	
	ATTRIBUTES:
		nextResDirection: Direction -- This is the relative direction of 
			the next residue in the chain.  If this residue is the C-terminus
			it has no next residue and the next direction is END of Direction
			enumerated type in Param.
		type: ResidueType -- This is the type of residue which can be any of
			the ResidueType enumerated class.
		terminus: TerminalType -- An additional residue type specifying if
			the residue is at the end of a chain, if so which end, or in the middle.
			Allows for end of chain residues to maintain a amino acid type a well
	CONSTRUCTORS:
		Residue(): No-argument -- Instantiates a residue with no next residue direction
			NULL type, and Not a terminus
		Residue(ResidueType type, TerminalType terminus): end residue type-specified 
			Constructor : sets the input ResidueType and TerminalType sets the 
			nextResDirection to END.
		Residue(Direction direction, ResidueType type, TerminalType terminus):
			direction and type specified Constructor -- sets ResidueType, TerminalType,
			and nextResDirection according to input parameters.
		
	METHODS:
	GETTERS
	SETTERS
	
	
	
	
	Version 1.0 beta 
		
  ------------------------------------------------------------------------------
 */


package protein_folder;

import protein_folder.Param.Direction;
import protein_folder.Param.ResidueType;
import protein_folder.Param.TerminalType;

public class Residue
{

//==================ATTRIBUTES======================================================
	private Direction nextResDirection;
	private ResidueType type;
	private TerminalType terminus;
//==================================================================================


	
//////__________Constructors_________________________________________________________

//----No-arg constructor used as the end residue-------------------------
public Residue(){
	this.nextResDirection = Param.END;
	this.type = Param.NULL;
	this.terminus = Param.NOT_TERMINUS;
}//end no arg constructor
//-----------------------------------------------------------------------

//----Argumented constructor input: (ResidueType, TerminalType)----------
public Residue(ResidueType inType, TerminalType inTerminus){
	this.terminus = inTerminus;
	this.nextResDirection = Param.END;
	this.type = inType;
}
//------------------------------------------------------------------------


//------Argumented Constructor input: (Direction,  ResidueType, TerminalType)-------
public Residue(Direction direction, ResidueType inType, TerminalType inTerminus){
	this.terminus = inTerminus;
	this.nextResDirection = direction;
	this.type = inType;
}
//----------------------------------------------------------------------
//___________________________________________________________________________________


//////++++++++GETTERS & SETTERS+++++++++++++++++++++++++++++++++++++++++++++++++++++

public ResidueType getType(){
	return type;
}

//get the direction of the residue next in the chain
public Direction getNextResDirection(){
	return this.nextResDirection;
};

//set the direction the next residue in the chain will be
public void setNextResDirection(Direction newDirection){
	nextResDirection = newDirection;
}

//get which terminal, if at all, the residue is
public TerminalType getTerminus() {
	return terminus;
}

//set the residue to be a terminus or not one
public void setTerminus(TerminalType terminus) {
this.terminus = terminus;
}
//////++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



}//end residue class