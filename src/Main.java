//import java.awt.Dimension;
import edu.uci.ics.jung.algorithms.scoring.ClosenessCentrality;
import edu.uci.ics.jung.algorithms.scoring.EigenvectorCentrality;
import edu.uci.ics.jung.algorithms.importance.BetweennessCentrality;
import java.util.Arrays;
import edu.uci.ics.jung.algorithms.cluster.*;
import java.util.Collection;
import java.util.Map;
import java.io.*;
import edu.uci.ics.jung.graph.*;
//import edu.uci.ics.jung.visualization.BasicVisualizationServer;
//import edu.uci.ics.jung.algorithms.layout.*;
//import javax.swing.*;
import java.util.Set;
import java.util.HashMap;
import java.util.Iterator;


class Link{
	public int Node1;
	public int Node2;
	
	Link(int inNode1, int inNode2){
		Node1 = inNode1;
		Node2 = inNode2;
	}
}

class ValueDegree{
	public double Degree;
	public int Element;
	
	ValueDegree(double x, int y){
		Degree=x;
		Element=y;
	}
}



public class Main {
	
	//Visualization
	/*public static void Visualize(Graph<Integer,Link> graph) {
		FRLayout2<Integer, Integer> layout = new FRLayout2(graph);
		layout.setSize(new Dimension(700,700)); // sets the initial size of the space
		BasicVisualizationServer<Integer,Integer> vv =
		new BasicVisualizationServer<Integer,Integer>(layout);
		vv.setPreferredSize(new Dimension(1920,1080)); //Sets the viewing area size

		JFrame frame = new JFrame("Simple Graph View");
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.getContentPane().add(vv);
		frame.pack();
		frame.setVisible(true); 
	}*/
	//Making Graph from TextFile
	public static void ReadTextFile(String fileName, Graph<Integer, Link> graph) {
		String line = null;
		BufferedReader buffRead = null;
		try {
			FileReader fRead = new FileReader(fileName);
			
			buffRead = new BufferedReader(fRead);
			
			while((line = buffRead.readLine()) != null) {
				String[] Nodes = line.split("\t");
				int node1 = Integer.parseInt(Nodes[0]);
				int node2 = Integer.parseInt(Nodes[1]);
				graph.addVertex(node1); graph.addVertex(node2);
				graph.addEdge(new Link(node1, node2), node1, node2);
			}
		} catch (FileNotFoundException ex) {
			System.out.println("File not found!");
		} catch (IOException ex) {
			System.out.println("Failed to read file!");
		} finally {
			try {
				if(buffRead != null) {
					buffRead.close();
				}
			} catch (IOException ex) {
				System.out.println("Failed to close BufferedReader!");
			}
		}
	}
	//Making Giant Connected Component 
	public static void MakeGCC(Graph<Integer, Link> graph) {
		WeakComponentClusterer<Integer, Link> wcc = new WeakComponentClusterer<Integer, Link>();
		Set<Set<Integer>> comps=wcc.transform(graph);
		Iterator<Set<Integer>>  cit = comps.iterator();
		int max=0;
		while (cit.hasNext()) {
			Set<Integer> compls= cit.next();
			if (max<=compls.size()) {
				max=compls.size();
			}
		}
		
		Iterator<Set<Integer>> iter = comps.iterator();
		while (iter.hasNext()) {
			Set<Integer> compls=iter.next();
			if (max!=compls.size()) {
				Iterator<Integer> its = compls.iterator();
				while (its.hasNext()) {
					Integer i=its.next();
					graph.removeVertex(i);
				}
			}
		}
	}
	//Pearson's Correlation Coefficient 
	public static double PearsonCorrelationCoefficient(double[] X,double[] Y) {
		double sumX=0;
		double sumY=0;
		int cnt=X.length;
		double[] X1=new double[cnt];
		double[] Y1=new double[cnt];
		for (int i=0; i<cnt;i++) {
			sumX=X[i]+sumX;
			sumY=Y[i]+sumY; 
			X1[i]=X[i];
			Y1[i]=Y[i];
		}
		double sumX1 = sumX;
		double cnt1 =cnt;
		double avX = sumX1/cnt1;
		double sumY1 = sumY;
		double avY = sumY1/cnt1;
		double coeff=0.0;
		double a=0.0;
		double b=0.0;
		double c=0.0;
		for (int i=0; i<cnt;i++) {
			a=a+(X1[i]-avX)*(Y1[i]-avY);
			b=b+(X1[i]-avX)*(X1[i]-avX);
			c=c+(Y1[i]-avY)*(Y1[i]-avY);
		}
		coeff=a/(Math.sqrt(b)*Math.sqrt(c));
	
	return coeff;
	}
	public static double Pearson(Graph<Integer, Link> graph) {

		
		double[] X = new double[graph.getEdgeCount()*2];
		double[] Y = new double[graph.getEdgeCount()*2];
		
		int cnt=0;
		 for (Link link: graph.getEdges()) {
	            int n1=link.Node1;
	            int n2=link.Node2;
	            
	            X[cnt]=graph.degree(n1);
	            Y[cnt]=graph.degree(n2);
	            ++cnt;
	            
	            X[cnt]=graph.degree(n2);
	            Y[cnt]=graph.degree(n1);
	            ++cnt;
		 }
	return PearsonCorrelationCoefficient(X,Y);
	}
	//Spearman's Correlation Coefficient
	public static double[] sort(double[] arr) {
		double pom=arr[0];
		for (int i=0; i<arr.length-1;i++) { 
            for (int j=i+1;j<arr.length;j++) {
              if (arr[i]>=arr[j]) {
                                pom=arr[i]; 
                                arr[i]=arr[j]; 
                                arr[j]=pom; 
              }
            }
		}
		for (int j=0;j<arr.length;j++) {
            if (arr[arr.length-1]<=arr[j]) {
                              pom=arr[arr.length-1]; 
                              arr[arr.length-1]=arr[j]; 
                              arr[j]=pom; 
            }
		}
	return arr;
	}
	public static double[] rang(double[] arr) {
		int br=0;
		int num=0;
		int sum=0;
		int j=0;
		double[] pom=arr;
		Arrays.sort(pom);
		double[] result=new double[pom.length];
		for(int i=0;i<pom.length-1;i++) {
			br+=1;
			if (pom[i]!=pom[i+1]) {
				result[i]=br;
			} else {
				num=br;
				sum=br;
				j=i;
				while ((j<pom.length-2) && (pom[j]==pom[j+1])) {
					num+=1;
					sum+=num;
					j+=1;
				}
				
				for(int k=i;k<num;k++) {
					double sum1=sum;
					double num1=num-1;
					result[k]=sum1/num1;
				}
				i=num-1;
				br=num;
			}
		
		}
		result[arr.length-1]=br+1;
		/*for (int i=0;i<result.length;i++) {
			System.out.println();
			System.out.print(result[i]);
			System.out.print(" ");
			System.out.print(arr[i]);
		}*/
	return result;
	}
	public static double Spearman(Graph<Integer, Link> graph) {
		double[] X = new double[graph.getEdgeCount()*2];
		double[] Y = new double[graph.getEdgeCount()*2];
		
		int cnt=0;
		 for (Link link: graph.getEdges()) {
	            int n1=link.Node1;
	            int n2=link.Node2;
	            
	            X[cnt]=graph.degree(n1);
	            Y[cnt]=graph.degree(n2);
	            ++cnt;
	            
	            X[cnt]=graph.degree(n2);
	            Y[cnt]=graph.degree(n1);
	            ++cnt;
		 }
		 
		double[] X1=new double[X.length];
		X1=rang(X);
		double[] Y1=new double[Y.length];
		Y1=rang(Y);
		return PearsonCorrelationCoefficient(X1,Y1);
	}
	/////////////////////////////GCC Dissolving///////////////////////////////////////////
	//Degree Dissolving
	public static void dissolveDegree(Graph<Integer, Link> graph) {
		ValueDegree[] X=new ValueDegree[graph.getVertexCount()];
		int k=0;
		Iterator<Integer> iter = graph.getVertices().iterator();
		while (iter.hasNext()) {
			Integer n = iter.next();
			X[k].Degree=graph.degree(n);
			X[k].Element=n;
			k++;
		}
		
		ValueDegree pom=X[0];
		for (int i=0; i<X.length-1;i++) { 
            for (int j=i+1;j<X.length;j++) {
              if (X[i].Degree>=X[j].Degree) {
                                pom=X[i]; 
                                X[i]=X[j]; 
                                X[j]=pom; 
              }
            }
		}
		
		for (int j=0;j>X.length;j++) {
            if (X[X.length-1].Degree>X[j].Degree) {
                              pom=X[X.length-1]; 
                              X[X.length-1]=X[j]; 
                              X[j]=pom; 
            }
		}
		
		for(int j=0; j<graph.getVertexCount();j++) {
			graph.removeVertex(X[j].Element);
		}
	}
	//BetweennessCentrality Dissolving
	public static void dissolveBetweennessCentrality(Graph<Integer, Link> graph) {
		
		ValueDegree[] X=new ValueDegree[graph.getVertexCount()];
		int k=0;
		BetweennessCentrality<Integer,Link> bc = new BetweennessCentrality<Integer,Link>(graph);
		Iterator<Integer> it=graph.getVertices().iterator();
		while (it.hasNext()) {
			Integer n=it.next();
			X[k].Element=n;
			X[k].Degree=bc.getVertexRankScore(n);
			k++;
		}
		
		ValueDegree pom=X[0];
		for (int i=0; i<X.length-1;i++) { 
            for (int j=i+1;j<X.length;j++) {
              if (X[i].Degree>=X[j].Degree) {
                                pom=X[i]; 
                                X[i]=X[j]; 
                                X[j]=pom; 
              }
            }
		}
		
		for (int j=0;j>X.length;j++) {
            if (X[X.length-1].Degree>X[j].Degree) {
                              pom=X[X.length-1]; 
                              X[X.length-1]=X[j]; 
                              X[j]=pom; 
            }
		}
		
		for(int j=0; j<graph.getVertexCount();j++) {
			graph.removeVertex(X[j].Element);
		}
		
	}
	//ClosenessCentrality Dissolving
	public static void dissolveClosenessCentrality(Graph<Integer,Link> graph) {
		ValueDegree[] X=new ValueDegree[graph.getVertexCount()];
		int k=0;
		ClosenessCentrality<Integer,Link> cc = new ClosenessCentrality<Integer,Link>(graph);
		Iterator<Integer> it=graph.getVertices().iterator();
		while (it.hasNext()) {
			Integer n=it.next();
			X[k].Element=n;
			X[k].Degree=cc.getVertexScore(n);
			k++;
		}
		
		ValueDegree pom=X[0];
		for (int i=0; i<X.length-1;i++) { 
            for (int j=i+1;j<X.length;j++) {
              if (X[i].Degree>=X[j].Degree) {
                                pom=X[i]; 
                                X[i]=X[j]; 
                                X[j]=pom; 
              }
            }
		}
		
		for (int j=0;j>X.length;j++) {
            if (X[X.length-1].Degree>X[j].Degree) {
                              pom=X[X.length-1]; 
                              X[X.length-1]=X[j]; 
                              X[j]=pom; 
            }
		}
		
		for(int j=0; j<graph.getVertexCount();j++) {
			graph.removeVertex(X[j].Element);
		}
		
	}
	//EigenvectorCentrality Dissolving
	public static void dissolveEigenvectorCentrality(Graph<Integer,Link> graph) {
		ValueDegree[] X=new ValueDegree[graph.getVertexCount()];
		int k=0;
		EigenvectorCentrality<Integer,Link> ec = new EigenvectorCentrality<Integer, Link>(graph);
		Iterator<Integer> it=graph.getVertices().iterator();
		while (it.hasNext()) {
			Integer n=it.next();
			X[k].Element=n;
			X[k].Degree=ec.getVertexScore(n);
			k++;
		}
		
		ValueDegree pom=X[0];
		for (int i=0; i<X.length-1;i++) { 
            for (int j=i+1;j<X.length;j++) {
              if (X[i].Degree>=X[j].Degree) {
                                pom=X[i]; 
                                X[i]=X[j]; 
                                X[j]=pom; 
              }
            }
		}
		
		for (int j=0;j>X.length;j++) {
            if (X[X.length-1].Degree>X[j].Degree) {
                              pom=X[X.length-1]; 
                              X[X.length-1]=X[j]; 
                              X[j]=pom; 
            }
		}
		
		for(int j=0; j<graph.getVertexCount();j++) {
			graph.removeVertex(X[j].Element);
		}
		
	}
	public static Graph<Integer, Link> copy(Graph<Integer,Link> graph) {
		Graph<Integer, Link> graphCopy = new UndirectedSparseGraph<Integer,Link>();
		for (Integer n: graph.getVertices()) {
			graphCopy.addVertex(n);
		}
		for (Link l: graph.getEdges()) {
			graphCopy.addEdge(l, graph.getIncidentVertices(l));
		}
	return graphCopy;
	}
	public static void main(String args[]) {
		
		Graph<Integer, Link> graph = new UndirectedSparseGraph<Integer, Link>();
		ReadTextFile("CA-GrQc.txt", graph);
		//Visualize(graph);
		MakeGCC(graph);
		Graph<Integer, Link> graphCopy = new UndirectedSparseGraph<Integer,Link>();
		graphCopy=copy(graph);
		
		System.out.println(Pearson(graph));
		System.out.println(Spearman(graph));
		//double[] k=new double[]{65, 32, 28,32,54,97,71};
		//k=rang(k);

	}
}
