package protein_folder;

public class PNG_Footer {

	final private String NI = "not init";
	private String solutionStr;
	private String energyStr;
	private String sequenceInfoStr;
	private String colInfoStr;
	
	public PNG_Footer(){
		solutionStr = NI;
		energyStr = NI;
		sequenceInfoStr = NI;
		colInfoStr = NI;
	}
	
	public String getSolutionStr() {
		return solutionStr;
	}
	public void setSolutionStr(String solutionStr) {
		this.solutionStr = solutionStr;
	}
	public String getEnergyStr() {
		return energyStr;
	}
	public void setEnergyStr(String energyStr) {
		this.energyStr = energyStr;
	}
	public String getSequenceInfoStr() {
		return sequenceInfoStr;
	}
	public void setSequenceInfoStr(String sequenceInfoStr) {
		this.sequenceInfoStr = sequenceInfoStr;
	}
	public String getColInfoStr() {
		return colInfoStr;
	}
	public void setColInfoStr(String colInfoStr) {
		this.colInfoStr = colInfoStr;
	}
	
}
