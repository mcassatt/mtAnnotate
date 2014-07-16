

public class genbankFeature implements Comparable<genbankFeature>{

		//contains type of feature
		public String type;
		//contains additional information about the feature
		//each datum is separated by \n.
		public String description;
		//start and end points of feature in the FASTA sequence
		//to be annotated
		public int[] location;
		
		public genbankFeature(String type, String des, int[] location){
			this.type = type;
			this.description =des;
			this.location = location;
		}
		//ascending sort by default, if same start points, then
		//the one with the first end point is less
		public int compareTo(genbankFeature feature){
				return this.location[0] - feature.location[0];
		
	}
}
