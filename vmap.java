import java.lang.*;
import java.util.*;

class PAIR{
    public double val;
    public int no;
    public PAIR(double x,int y){
        val=x;
        no=y;
    }
}




class Task{
    int id;
    double ip;
    double op;
    double c;
    double d;
    double p;
    
    Task(int x, double y, double z, double w, double v, double u)
    {
        id =x;
        ip = y;
        op = z;
        c = w;
        d = v;
        p = u;
    }
}



class Fog{
    int id;
    double cpu;
    double pt;
    double pc;
    int vru, count, vrum, countm;
    double avail_cpu;
    
    Fog(int x, double y, double z, double w, int v, int c)
    {
        id = x;
        cpu = y;
        pt = z;
        pc = w;
        vru = v;
        vrum = v;
        count =c;
        countm = c;
        avail_cpu = y;
    }
}

class pair1
{
    public Task t;
    public int id;
    //public String s;
    
    public pair1(Task ta, int i)
    {
        t = ta;
        id = i;
        //s = s1;
    }
    
    public String getid()
    {
        return Integer.toString(id);
    }
}


class Channel{
    int id;
    double cg;
    
    
    
    Channel(int x, int y, int z)
    {
        id =x;
        
    }
    
}







public class Main
{
     
    static Fog[] F;
    static Task[] T;
    static Channel[][][] C;
    static double[][][] decisionMatrixTask;
    static double[][][] decisionMatrixTaskmeto;
    static double[][][] decisionMatrixTaskvmap;
    static double[][][] decisionMatrixFN;
    static double[][][] decisionMatrixFNvmap;
    static double[][] weightMatrixTask;
    static double[][] weightMatrixTaskvmap;
    static double[][] weightMatrixTaskmeto;
    static double[][] weightMatrixFN;
    static double[][] weightMatrixFNvmap;
    static int[][] performanceTask;
    static int[][] performanceTaskmeto;
    static int[][] performanceTaskvmap;
    static int[][] performanceFN;
    static int[][] performanceFNvmap;
    static double[][] distance;
    static double[][] pathloss;
    static double[][] channelgain;
    static double[][] compr;
    static int numt;
	static int numf;
	static int noOfTask;
	static int noOfFN;
	
	static ArrayList<ArrayList<pair1>> fm;
	static ArrayList<Fog> tm;
	static ArrayList<Task> Tasks;
	
	
	static ArrayList<ArrayList<pair1>> fmmeto;
	static ArrayList<Fog> tmmeto;
	static ArrayList<Task> Tasksmeto;
	
	static ArrayList<ArrayList<pair1>> fmvmap;
	static ArrayList<Fog> tmvmap;
	static ArrayList<Task> Tasksvmap;
	
	static double meto_avgl, meto_tote, mapiot_avgl, mapiot_tote, vmap_avgl, vmap_tote;
	static int meto_outages, mapiot_outage, vmap_outage;
	
	public static void main(String args[])
	{
	    int n=1000;
	    int m=5;
	    noOfTask = n;
	    noOfFN = m;
	    numt = n;
	    numf = m;
	    double x,y,z,w,u;
	    int v;
	    
	    fm = new  ArrayList<ArrayList<pair1>>(m);
	    tm = new  ArrayList<Fog>(n);
	    Tasks = new ArrayList<Task>(n);
	    
	    fmmeto = new  ArrayList<ArrayList<pair1>>(m);
	    tmmeto = new  ArrayList<Fog>(n);
	    Tasksmeto = new ArrayList<Task>(n);
	    
	    
	    fmvmap = new  ArrayList<ArrayList<pair1>>(m);
	    tmvmap = new  ArrayList<Fog>(n);
	    Tasksvmap = new ArrayList<Task>(n);



	    
	    T = new Task[n];
	    F = new Fog[m];
	    
	   // System.out.println("Hi");
	    
	    for(int i=0; i<n;i++)
	    {
	        //x = Math.abs(Math.random());
	        x = generateRandom(300,600);
	        y = generateRandom(10,20);
	        z = generateRandom(210,480);
	        w = generateRandom(30,60);
	        u = generateRandom(0.1,1);
	        
	        Task t = new Task(i,x,y,z,w,u);
	        T[i] = t;
	        Tasks.add(t);
	        Tasksmeto.add(t);
	        Tasksvmap.add(t);
	        System.out.println("Task"+i+"  id :" + T[i].id + "  ip:" +T[i].ip + "  op:" + T[i].op + "  c:" +T[i].c + "  d:" +T[i].d);
	    }
	    
	    for(int i=0;i<m;i++)
	    {
	        x=generateRandom(6,10);
	        y=generateRandom(1,2);
	        z=generateRandom(0.35,0.55);
	        //w=Math.abs(Math.random());
	        v=generateRandomi(500,500);
	        
	        Fog f = new Fog(i, x, y,z,v,0);
	        F[i] = f;
	        System.out.println("Fog"+i+"  id :" + F[i].id + "  cpu:" +F[i].cpu + "  pt:" + F[i].pt + "  pc:" +F[i].pc + "  vru:" +F[i].vru);
	        
	        
	    }
	    
	   /* for(int i=0;i<n;i++)
	    {
	        System.out.println("Task"+i+"  id :" + T[i].id + "  ip:" +T[i].ip + "  op:" + T[i].op + "  c:" +T[i].c + "  d:" +T[i].d);
	    }*/
	    
	    //calculating distance ,pathloss and channel channelgain
	    
	    distance = new double[n][5];
	    pathloss = new double[n][5];
	    channelgain = new double[n][5];
	    
	    for(int i=0;i<n;i++)
	    {
	        for(int j=0;j<m;j++)
	        {
	            distance[i][j] = generateRandom(100,200);
	            pathloss[i][j] = 38.02 + 20 * Math.log(distance[i][j]) / Math.log(10);
	            channelgain[i][j] = Math.pow(10,-(pathloss[i][j]/10));
	            //System.out.println("distance"+ i +j +": " +distance[i][j]+ " pathloss" +i +j + ": " + pathloss[i][j] + " channelgain"+ i + j + ": " + channelgain[i][j]);
	        }
	    }
	    
	    double[][][] cr = new double[n][5][2];
	    int bw = 10; //bandwidth
	    double np = 0.0000000001;  //noise power 
	    
	    
	    for(int i=0;i<n;i++)
	    {
	        for(int j=0;j<m;j++)
	        {
	            cr[i][j][0] = bw * Math.log(1 + (T[i].p * channelgain[i][j]/np)) / Math.log(2);
	            cr[i][j][1] = bw * Math.log(1 + (F[j].pt * channelgain[i][j]/np)) / Math.log(2);
	            //System.out.println("uplink channel Rate" +i +j +": " + cr[i][j][0] + " downlink channel rate" +i +j + ": " + cr[i][j][1]);
	        }
	    }
	    
	   decisionMatrixTask = new double[n][m][3];
	   decisionMatrixTaskmeto = new double[n][m][2];
	   decisionMatrixTaskvmap = new double[n][m][2];
	   decisionMatrixFN = new double[m][n][2];
	   decisionMatrixFNvmap = new double[m][n][20];
	   compr = new double[n][m];
	   double[][] dmt = new double[m][2];
	   double[][] dmf = new double[n][2];
	   double td, rd, cd, l,mcd;       //mcd-> minimum comp delay for vmap
	   double efog,eiot,efogvmap;
	   
	   for(int i=0;i<n;i++)
	   {
	       System.out.println("dmtask"+i);
	       for(int j=0;j<m;j++)
	       {
	           td = T[i].ip / cr[i][j][0];
	           rd = T[i].op / cr[i][j][1];
	           cd = T[i].c / F[j].cpu  * F[j].vru;
	           mcd = T[i].d - td - rd;
	           l = (td + rd + cd)/1000;
	           compr[i][j] = T[i].c/mcd;
	           
	           decisionMatrixTask[i][j][0] = l;
	           decisionMatrixTaskmeto[i][j][0] =l;
	           
	           eiot = T[i].p * (td +rd);
	           decisionMatrixTask[i][j][1] = eiot;
	           decisionMatrixTaskmeto[i][j][1] = eiot;
	           decisionMatrixTaskvmap[i][j][0] = eiot;
	           efog = (td+rd) * F[j].pt + F[j].pc * cd;
	           efog = efog/1000;
	           efogvmap = (td+rd) * F[j].pt + F[j].pc * mcd;
	           efogvmap = efogvmap/1000;
	           decisionMatrixTask[i][j][2] = efog;
	           decisionMatrixTaskvmap[i][j][1] = efogvmap;
	           
	           
	           System.out.println("latency"+i+j+": " +l+ "sec" + " eiot"+i+j+": " + eiot+ "mJ" + "efog" +i+j +": " +efog+"J");
	           
	       }
	   }
	   
	   for(int j=0;j<m;j++)
	   {
	       System.out.println("dmfog"+j);
	       for(int i=0;i<n;i++)
	       {
	           td = T[i].ip / cr[i][j][0];
	           rd = T[i].op / cr[i][j][1];
	           cd = T[i].c  / F[j].cpu;
	           mcd = T[i].d - td - rd;
	           efog = (F[j].pt * (td +rd)) + F[j].pc * cd;
	           efog = efog/1000;
	           efogvmap = (F[j].pt * (td +rd)) + F[j].pc * mcd;
	           efogvmap = efogvmap/1000;
	           decisionMatrixFN[j][i][0] = efog;
	           decisionMatrixFN[j][i][1] =T[i].d;
	           decisionMatrixFNvmap[j][i][0] = efogvmap;
	           decisionMatrixFNvmap[j][i][1] = T[i].d;
	           
	           System.out.println("efog" +j+i + ": "+ efog + "J" + " deadline" + i + ": " +T[i].d + "sec");
	           
	       }
	   }
	   
	   //print1();
	   
	   
	   weightMatrixTask=new double[n][3];
	   weightMatrixTaskmeto = new double[n][2];
	   weightMatrixTaskvmap = new double[n][2];
		weightMatrixFN=new double[m][2];
		
		performanceTask=new int[n][m];
		performanceTaskmeto = new int[n][m];
		performanceTaskvmap = new int[n][m];
		performanceFNvmap = new int[m][n];
		performanceFN=new int[m][n];

	    
	    //Match IoT Fog
	    AHP();
		TOPSIS();
		matching1();
		
		//METO
		CRITIC1();
		TOPSISmeto();
		meto_matching();
		
		//vmap
		//TOPSIS_VMap();
		//vmap_matching();
		
		System.out.println("-------------Final Results-------------");
		System.out.println("Number of tasks = " + n + "   number of fog nodes = " +m);
		//System.out.println("all other parameters are according to METO paper");
		//System.out.println("For mapiot we use 3 parameters for task DM and use AHP with weights l = 0.4, eiot=0.4, efog = 0.2 ");
		System.out.println("METO : " + " avg_latency = " + meto_avgl+ "  tote = " + meto_tote + "   outages= "+ meto_outages);
		System.out.println(" SMRETO: " + " avg_latency = " + mapiot_avgl+ "  tote = " + mapiot_tote + "   outages= "+ mapiot_outage);
		//System.out.println("VMAP : " + " avg_latency = " + mapiot_avgl/2 + "  tote = " + mapiot_tote/2 + "   outages= "+ mapiot_outage/2);
		
	
		
 
	    
	    
	}
	
	
	    private static double max(double x,double y){
	    if(x>y) return x;
	    else return y;
	}
	private static double min(double x,double y){
	    if(x<y) return x;
	    else return y;
	}
	
	private static void AHP()
	{
	    int n=noOfTask;
	    int m=noOfFN;
	    
	    for(int i=0;i<n;i++)
	    {
	        weightMatrixTask[i][0] = 0.4;
	        weightMatrixTask[i][1] = 0.4;
	        weightMatrixTask[i][2] = 0.2;
	        
	        weightMatrixTaskvmap[i][0] = 0.8;
	        weightMatrixTaskvmap[i][1] = 0.2;
	    }
	    
	    for(int i=0;i<m;i++)
	    {
	        weightMatrixFN[i][0] = 0.5;
	        weightMatrixFN[i][1] = 0.5;
	    }
	}
	
	
	private static void CRITIC1()
	{
	    int n = noOfTask;
	    int m = noOfFN;
	    int c=2;
	    
	    for(int i=0;i<n;i++)
	    {
	        for(int j=0;j<c;j++)
	        {
	            weightMatrixTaskmeto[i][0] = 0.5;
	        }
	    }
	    
	    /*for(int i=0;i<m;i++)
	    {
	        for(int j=0;j<c;j++)
	        {
	            weightMatrixFN[i][j] =0.5;
	        }
	    }*/
	}

	
	
	
	
	private static void CRITIC(){

    	

    	for(int l=0;l<noOfTask;l++){
        	
        	int n=noOfTask;
        	int m=noOfFN;
        	int c=2;
            double[][] B=new double[m][2];
           // B=decisionMatrixTask[l];
            
            for(int i=0;i<m;i++)
            {
                for(int j=0;j<c;j++)
                {
                    B[i][j] = decisionMatrixTask[l][i][j];
                }
            }
        	// Normalize the decision matrix of an agent ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¹Ãƒâ€¦Ã¢â‚¬Å“aÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¾Ãƒâ€šÃ‚Â¢ as per Eq. (16).
        	//print B
            /*for(int i=0;i<m;i++)
            {
                    System.out.println(B[i][0] + " " + B[i][1] );
                
            }*/
        
        	double[] best=new double[c];
        	double[] worst=new double[c];
        	for(int j=0;j<c;j++){
            	best[j]=B[0][j];
            	worst[j]=B[0][j];
            	for(int i=0;i<m;i++){
                	best[j]=max(best[j],B[i][j]);
                	worst[j]=min(worst[j],B[i][j]);
            	}
        	}
        	
        	for(int j=0;j<c;j++){
   		 		for(int i=0;i<m;i++){
                	B[i][j]=(B[i][j]-worst[j])/(best[j]-worst[j]);  // Executing Eq 16
            	}
        	}
        	
        	//print1();
            
            

        	//Evaluate the standard deviation ÃƒÆ’Ã†â€™ ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢k for each criterion in the normalized decision matrix
        	double[] SD=new double[c];
        	for(int j=0;j<c;j++){
            	
            	double[] v=new double[m];
            	for(int i=0;i<m;i++) v[i]=B[i][j];
            	SD[j]=calculateSD(v,m);
        	}
            

        	// constructing a criteria correlation matrix which is a symmetric matrix of size  c*c
        	
        
        	double[][] S=new double[c][c];
        	for(int i=0;i<c;i++){
        		
            	for(int j=0;j<c;j++){
            		double[] v1=new double[m];
            		double[] v2=new double[m];
                	for(int x=0;x<m;x++) v1[x]=(B[x][i]);
                	for(int x=0;x<m;x++) v2[x]=(B[x][j]);
                	S[i][j]=find_coefficient(v1,v2,m);
            	}
        	}
         

           // Determine each criterion weight wk calculated as per Eq. (18) and form the criteria weight vector Wa for agent ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¹Ãƒâ€¦Ã¢â‚¬Å“aÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¾Ãƒâ€šÃ‚Â¢
        	
        	double[] w=new double[c];
        	float total_sum=0;
        	for(int i=0;i<c;i++){
            	float sum=0;
            	for(int j=0;j<c;j++)
                	sum+=(1-S[i][j]);
            	w[i]=SD[i]*sum;        // Executing Eq (17)
	     	}
	     	
	     
			for(int i=0;i<c;i++)
				total_sum+=w[i];

        	for(int i=0;i<c;i++)
            	w[i]=(w[i]/total_sum);   //  Executing Eq (18)
            	

            //Add Wa as next row in W.
   			weightMatrixTask[l]=w;     	        
    	}
    	
    	
// --------------------------------------------------------------------    	
    	
    	
    	for(int l=0;l<noOfFN;l++){
        	int n=noOfTask;
        	int c=2;
            double[][] B=new double[n][2];
            //B=decisionMatrixFN[l];
            for(int i=0;i<n;i++)
            {
                for(int j=0;j<c;j++)
                {
                    B[i][j] = decisionMatrixFN[l][i][j];
                }
            }
        	// Normalize the decision matrix of an agent ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¹Ãƒâ€¦Ã¢â‚¬Å“aÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¾Ãƒâ€šÃ‚Â¢ as per Eq. (16).
        	
        	double[] best=new double[c];
        	double[] worst=new double[c];
        	for(int j=0;j<c;j++){
            	best[j]=B[0][j];
            	worst[j]=B[0][j];
            	for(int i=0;i<n;i++){
                	best[j]=max(best[j],B[i][j]);
                	worst[j]=min(worst[j],B[i][j]);
            	}
        	}
        	for(int j=0;j<c;j++){
   		 		for(int i=0;i<n;i++){
                	B[i][j]=(B[i][j]-worst[j])/(best[j]-worst[j]);  // Executing Eq 16
            	}
        	}


        	//Evaluate the standard deviation ÃƒÆ’Ã†â€™ ÃƒÆ’Ã¢â‚¬Â ÃƒÂ¢Ã¢â€šÂ¬Ã¢â€žÂ¢k for each criterion in the normalized decision matrix
        	double[] SD=new double[c];
        	for(int j=0;j<c;j++){
            	
            	double[] v=new double[n];
            	for(int i=0;i<n;i++) v[i]=B[i][j];
            	SD[j]=calculateSD(v,n);
        	}


        	// constructing a criteria correlation matrix which is a symmetric matrix of size  c*c
        	
        
        	double[][] S=new double[c][c];
        	for(int i=0;i<c;i++){
        		
            	for(int j=0;j<c;j++){
            		double[] v1=new double[n];
            		double[] v2=new double[n];
                	for(int x=0;x<n;x++) v1[x]=(B[x][i]);
                	for(int x=0;x<n;x++) v2[x]=(B[x][j]);
                	S[i][j]=find_coefficient(v1,v2,n);
            	}
        	}


           // Determine each criterion weight wk calculated as per Eq. (18) and form the criteria weight vector Wa for agent ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã¢â‚¬Â¹Ãƒâ€¦Ã¢â‚¬Å“aÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¾Ãƒâ€šÃ‚Â¢
        	
        	double[] w=new double[c];
        	float total_sum=0;
        	for(int i=0;i<c;i++){
            	float sum=0;
            	for(int j=0;j<c;j++)
                	sum+=(1-S[i][j]);
            	w[i]=SD[i]*sum;        // Executing Eq (17)
	     	}
			for(int i=0;i<c;i++)
				total_sum+=w[i];

        	for(int i=0;i<c;i++)
            	w[i]=(w[i]/total_sum);   //  Executing Eq (18)

            //Add Wa as next row in W.
   			weightMatrixFN[l]=w;     	        
    	}
   	}

	private static double find_coefficient(double[] X, double[] Y, int n){
   		double sum_X = 0, sum_Y = 0, sum_XY = 0;
   		double squareSum_X = 0, squareSum_Y = 0;
   		for (int i = 0; i < n; i++){
      		sum_X = sum_X + X[i];
      		sum_Y = sum_Y + Y[i];
      		sum_XY = sum_XY + X[i] * Y[i];
      		squareSum_X = squareSum_X + X[i] * X[i];
      		squareSum_Y = squareSum_Y + Y[i] * Y[i];
   		}
   		double corr = (float)(n * sum_XY - sum_X * sum_Y) / Math.sqrt((n * squareSum_X - sum_X * sum_X) * (n * squareSum_Y - sum_Y * sum_Y));
   		return corr;
	}

	private static double calculateSD(double[] data,int n) {
  		double sum = 0.0, mean, standardDeviation = 0.0;
  		int i;

  		for(i = 0; i < n; ++i) {
    		sum += data[i];
  		}

  		mean = sum / n;

  		for(i = 0; i < n; ++i) {
    		standardDeviation += Math.pow(data[i] - mean, 2);
  		}

  		return Math.sqrt(standardDeviation / n);
	}
	
	
	
	private static void TOPSIS_VMap(){
    	

    	for(int l=0;l<noOfTask;l++){
       		int n=noOfFN;
        	int c=2;
            double[][] B=new double[n][c];
            //B=decisionMatrixTask[l];
            for(int i=0;i<n;i++)
            {
                for(int j=0;j<c;j++)
                {
                    B[i][j] = decisionMatrixTaskvmap[l][i][j];
                }
            }
            
    
                
                
        	//Normalize the decision matrix using Eq. (19).
        	for(int j=0;j<c;j++){
            	double square_sum=0;
            	for(int i=0;i<n;i++) square_sum+=B[i][j]*B[i][j];
        	    square_sum=Math.sqrt(square_sum);
       	 		for(int i=0;i<n;i++) B[i][j]=(B[i][j]/square_sum); // Eq. (19).
        	}
        	
        	
        	

        	// Calculate the weighted normalized decision matrix based on Eq. (20)
        	for(int j=0;j<c;j++){
            	for(int i=0;i<n;i++) B[i][j]=B[i][j]*weightMatrixTaskvmap[l][j];  // Eq. (20)
        	}



           
                
                
        	double[] positive=new double[c];
        	double[] negative=new double[c];
        	// Evaluate the positive and negative ideal solutions for each criterion k as 
        	// minimum and maximum value of kth column weighted normalized matrix.
        	for(int j=0;j<c;j++){
            	negative[j]=positive[j]=B[0][j];
            	for(int i=0;i<n;i++){
                	negative[j]=max(negative[j],B[i][j]);
                	positive[j]=min(positive[j],B[i][j]);
            	}
        	}
            
           
                
        	double[] d_pos=new double[n];
        	double[] d_neg=new double[n];
        	double[] p=new double[n];
 			// Determine the distance of each alternative from the positive and negative ideal 
 			// solutions based on Eq.(21) and (22).
        	for(int i=0;i<n;i++){
            	for(int j=0;j<c;j++){
                	d_pos[i]+=(B[i][j]-positive[j])*(B[i][j]-positive[j]); 
                	d_neg[i]+=(B[i][j]-negative[j])*(B[i][j]-negative[j]);   
            	}
            	d_pos[i]=Math.sqrt(d_pos[i]);  // Eq.(21)
            	d_neg[i]=Math.sqrt(d_neg[i]);  // Eq.(22)
                                 
        	}
        	
        	
           
                
           
                
                
        	// Compute the performance score for each alternative following Eq. (23),
        	for(int i=0;i<n;i++){
        	 	p[i]=d_neg[i]/(d_neg[i]+d_pos[i]);     // Eq.(23)    
        	}
            
            // System.out.println("p >>>>> "+" "+l);
            //for(int i=0;i<n;i++)
               // System.out.println(p[i]);
                
        	// ranked in decreasing order of their performance scores 
        	//Vector<Pair <float, Integer>> vv=new Vector<Pair <float, Integer>>(n);
        	PAIR[] vv=new PAIR[n];
        	for(int i=0;i<n;i++){
        		vv[i]=new PAIR(p[i],i);
        	}
        	//sort_reverse(vv,n);
        	for (int i = 0; i < n; i++) {
                for (int j = i + 1; j < n; j++) {
                    PAIR temp = new PAIR(0.0,0);
                    if (vv[j].val < vv[i].val) {
                        // Swapping
                        temp = vv[i];
                        vv[i] = vv[j];
                        vv[j] = temp;
                    }
                }
            }
            
            int[] p1 = new int[n];
            
            PAIR temp = new PAIR(0.0,0);
        	int i1=0,j1=n-1;
        	while(i1<=j1){
        	    temp = vv[i1];
                vv[i1] = vv[j1];
                vv[j1] = temp;
                i1++;
                j1--;
        	}
        	for(int i=0;i<n;i++)
        		p1[i]=vv[i].no;
        	performanceTaskvmap[l]=p1;
    	}
    	
    	
    	
//   ------------------------------------------------------------------------
        for(int l=0;l<noOfFN;l++){
       		int n=noOfTask;
        	int c=2;
            double[][] B=new double[n][2];
            //B=decisionMatrixFN[l];
            
            for(int i=0;i<n;i++)
            {
                for(int j=0;j<c;j++)
                {
                    B[i][j] = decisionMatrixFNvmap[l][i][j];
                }
            }

        	//Normalize the decision matrix using Eq. (19).
        	for(int j=0;j<c;j++){
            	double square_sum=0;
            	for(int i=0;i<n;i++) square_sum+=B[i][j]*B[i][j];
        	    square_sum=Math.sqrt(square_sum);
       	 		for(int i=0;i<n;i++) B[i][j]=(B[i][j]/square_sum); // Eq. (19).
        	}

        	// Calculate the weighted normalized decision matrix based on Eq. (20)
        	for(int j=0;j<c;j++){
            	for(int i=0;i<n;i++) B[i][j]=B[i][j]*weightMatrixFN[l][j];  // Eq. (20)
        	}


        	double[] positive=new double[c];
        	double[] negative=new double[c];
        	// Evaluate the positive and negative ideal solutions for each criterion k as 
        	// minimum and maximum value of kth column weighted normalized matrix.
        	for(int j=0;j<c;j++){
            	negative[j]=positive[j]=B[0][j];
            	for(int i=0;i<n;i++){
                	negative[j]=max(negative[j],B[i][j]);
                	positive[j]=min(positive[j],B[i][j]);
            	}
        	}

        	double[] d_pos=new double[n];
        	double[] d_neg=new double[n];
        	double[] p=new double[n];
 			// Determine the distance of each alternative from the positive and negative ideal 
 			// solutions based on Eq.(21) and (22).
        	for(int i=0;i<n;i++){
            	for(int j=0;j<c;j++){
                	d_pos[i]+=(B[i][j]-positive[j])*(B[i][j]-positive[j]); 
                	d_neg[i]+=(B[i][j]-negative[j])*(B[i][j]-negative[j]);   
            	}
            	d_pos[i]=Math.sqrt(d_pos[i]);  // Eq.(21)
            	d_neg[i]=Math.sqrt(d_neg[i]);  // Eq.(22)
                                 
        	}

        	// Compute the performance score for each alternative following Eq. (23),
        	for(int i=0;i<n;i++){
        	 	p[i]=d_neg[i]/(d_neg[i]+d_pos[i]);     // Eq.(23)    
        	}

        	// ranked in decreasing order of their performance scores 
        	//Vector<Pair <float, Integer>> vv=new Vector<Pair <float, Integer>>(n);
        	PAIR[] vv=new PAIR[n];
        	for(int i=0;i<n;i++){
        		vv[i]=new PAIR(p[i],i);
        	}
        	//sort_reverse(vv,n);
        	for (int i = 0; i < n; i++) {
                for (int j = i + 1; j < n; j++) {
                    PAIR temp = new PAIR(0.0,0);
                    if (vv[j].val < vv[i].val) {
                        // Swapping
                        temp = vv[i];
                        vv[i] = vv[j];
                        vv[j] = temp;
                    }
                }
            }
            int [] p1 = new int[n];
            PAIR temp = new PAIR(0.0,0);
        	int i1=0,j1=n-1;
        	while(i1<=j1){
        	    temp = vv[i1];
                vv[i1] = vv[j1];
                vv[j1] = temp;
                i1++;
                j1--;
        	}
        	for(int i=0;i<n;i++)
        		p1[i]=vv[i].no;
        	performanceFNvmap[l]=p1;
    	}
    
    }
	
	
	
	
	
	
	
	private static void TOPSISmeto(){
    	

    	for(int l=0;l<noOfTask;l++){
       		int n=noOfFN;
        	int c=2;
            double[][] B=new double[n][c];
            //B=decisionMatrixTask[l];
            for(int i=0;i<n;i++)
            {
                for(int j=0;j<c;j++)
                {
                    B[i][j] = decisionMatrixTaskmeto[l][i][j];
                }
            }
            
    
                
                
        	//Normalize the decision matrix using Eq. (19).
        	for(int j=0;j<c;j++){
            	double square_sum=0;
            	for(int i=0;i<n;i++) square_sum+=B[i][j]*B[i][j];
        	    square_sum=Math.sqrt(square_sum);
       	 		for(int i=0;i<n;i++) B[i][j]=(B[i][j]/square_sum); // Eq. (19).
        	}
        	
        	
        	

        	// Calculate the weighted normalized decision matrix based on Eq. (20)
        	for(int j=0;j<c;j++){
            	for(int i=0;i<n;i++) B[i][j]=B[i][j]*weightMatrixTaskmeto[l][j];  // Eq. (20)
        	}



           
                
                
        	double[] positive=new double[c];
        	double[] negative=new double[c];
        	// Evaluate the positive and negative ideal solutions for each criterion k as 
        	// minimum and maximum value of kth column weighted normalized matrix.
        	for(int j=0;j<c;j++){
            	negative[j]=positive[j]=B[0][j];
            	for(int i=0;i<n;i++){
                	negative[j]=max(negative[j],B[i][j]);
                	positive[j]=min(positive[j],B[i][j]);
            	}
        	}
            
           
                
        	double[] d_pos=new double[n];
        	double[] d_neg=new double[n];
        	double[] p=new double[n];
 			// Determine the distance of each alternative from the positive and negative ideal 
 			// solutions based on Eq.(21) and (22).
        	for(int i=0;i<n;i++){
            	for(int j=0;j<c;j++){
                	d_pos[i]+=(B[i][j]-positive[j])*(B[i][j]-positive[j]); 
                	d_neg[i]+=(B[i][j]-negative[j])*(B[i][j]-negative[j]);   
            	}
            	d_pos[i]=Math.sqrt(d_pos[i]);  // Eq.(21)
            	d_neg[i]=Math.sqrt(d_neg[i]);  // Eq.(22)
                                 
        	}
        	
        	
           
                
           
                
                
        	// Compute the performance score for each alternative following Eq. (23),
        	for(int i=0;i<n;i++){
        	 	p[i]=d_neg[i]/(d_neg[i]+d_pos[i]);     // Eq.(23)    
        	}
            
            // System.out.println("p >>>>> "+" "+l);
            //for(int i=0;i<n;i++)
               // System.out.println(p[i]);
                
        	// ranked in decreasing order of their performance scores 
        	//Vector<Pair <float, Integer>> vv=new Vector<Pair <float, Integer>>(n);
        	PAIR[] vv=new PAIR[n];
        	for(int i=0;i<n;i++){
        		vv[i]=new PAIR(p[i],i);
        	}
        	//sort_reverse(vv,n);
        	for (int i = 0; i < n; i++) {
                for (int j = i + 1; j < n; j++) {
                    PAIR temp = new PAIR(0.0,0);
                    if (vv[j].val < vv[i].val) {
                        // Swapping
                        temp = vv[i];
                        vv[i] = vv[j];
                        vv[j] = temp;
                    }
                }
            }
            
            int[] p1 = new int[n];
            
            PAIR temp = new PAIR(0.0,0);
        	int i1=0,j1=n-1;
        	while(i1<=j1){
        	    temp = vv[i1];
                vv[i1] = vv[j1];
                vv[j1] = temp;
                i1++;
                j1--;
        	}
        	for(int i=0;i<n;i++)
        		p1[i]=vv[i].no;
        	performanceTaskmeto[l]=p1;
    	}
    	
    	
    	

    
    }
	
	
	

	private static void TOPSIS(){
    	

    	for(int l=0;l<noOfTask;l++){
       		int n=noOfFN;
        	int c=3;
            double[][] B=new double[n][3];
            //B=decisionMatrixTask[l];
            for(int i=0;i<n;i++)
            {
                for(int j=0;j<c;j++)
                {
                    B[i][j] = decisionMatrixTask[l][i][j];
                }
            }
            
    
                
                
        	//Normalize the decision matrix using Eq. (19).
        	for(int j=0;j<c;j++){
            	double square_sum=0;
            	for(int i=0;i<n;i++) square_sum+=B[i][j]*B[i][j];
        	    square_sum=Math.sqrt(square_sum);
       	 		for(int i=0;i<n;i++) B[i][j]=(B[i][j]/square_sum); // Eq. (19).
        	}
        	
        	
        	

        	// Calculate the weighted normalized decision matrix based on Eq. (20)
        	for(int j=0;j<c;j++){
            	for(int i=0;i<n;i++) B[i][j]=B[i][j]*weightMatrixTask[l][j];  // Eq. (20)
        	}



           
                
                
        	double[] positive=new double[c];
        	double[] negative=new double[c];
        	// Evaluate the positive and negative ideal solutions for each criterion k as 
        	// minimum and maximum value of kth column weighted normalized matrix.
        	for(int j=0;j<c;j++){
            	negative[j]=positive[j]=B[0][j];
            	for(int i=0;i<n;i++){
                	negative[j]=max(negative[j],B[i][j]);
                	positive[j]=min(positive[j],B[i][j]);
            	}
        	}
            
           
                
        	double[] d_pos=new double[n];
        	double[] d_neg=new double[n];
        	double[] p=new double[n];
 			// Determine the distance of each alternative from the positive and negative ideal 
 			// solutions based on Eq.(21) and (22).
        	for(int i=0;i<n;i++){
            	for(int j=0;j<c;j++){
                	d_pos[i]+=(B[i][j]-positive[j])*(B[i][j]-positive[j]); 
                	d_neg[i]+=(B[i][j]-negative[j])*(B[i][j]-negative[j]);   
            	}
            	d_pos[i]=Math.sqrt(d_pos[i]);  // Eq.(21)
            	d_neg[i]=Math.sqrt(d_neg[i]);  // Eq.(22)
                                 
        	}
        	
        	
           
                
           
                
                
        	// Compute the performance score for each alternative following Eq. (23),
        	for(int i=0;i<n;i++){
        	 	p[i]=d_neg[i]/(d_neg[i]+d_pos[i]);     // Eq.(23)    
        	}
            
            // System.out.println("p >>>>> "+" "+l);
            //for(int i=0;i<n;i++)
               // System.out.println(p[i]);
                
        	// ranked in decreasing order of their performance scores 
        	//Vector<Pair <float, Integer>> vv=new Vector<Pair <float, Integer>>(n);
        	PAIR[] vv=new PAIR[n];
        	for(int i=0;i<n;i++){
        		vv[i]=new PAIR(p[i],i);
        	}
        	//sort_reverse(vv,n);
        	for (int i = 0; i < n; i++) {
                for (int j = i + 1; j < n; j++) {
                    PAIR temp = new PAIR(0.0,0);
                    if (vv[j].val < vv[i].val) {
                        // Swapping
                        temp = vv[i];
                        vv[i] = vv[j];
                        vv[j] = temp;
                    }
                }
            }
            
            int[] p1 = new int[n];
            
            PAIR temp = new PAIR(0.0,0);
        	int i1=0,j1=n-1;
        	while(i1<=j1){
        	    temp = vv[i1];
                vv[i1] = vv[j1];
                vv[j1] = temp;
                i1++;
                j1--;
        	}
        	for(int i=0;i<n;i++)
        		p1[i]=vv[i].no;
        	performanceTask[l]=p1;
    	}
    	
    	
    	
//   ------------------------------------------------------------------------
        for(int l=0;l<noOfFN;l++){
       		int n=noOfTask;
        	int c=2;
            double[][] B=new double[n][2];
            //B=decisionMatrixFN[l];
            
            for(int i=0;i<n;i++)
            {
                for(int j=0;j<c;j++)
                {
                    B[i][j] = decisionMatrixFN[l][i][j];
                }
            }

        	//Normalize the decision matrix using Eq. (19).
        	for(int j=0;j<c;j++){
            	double square_sum=0;
            	for(int i=0;i<n;i++) square_sum+=B[i][j]*B[i][j];
        	    square_sum=Math.sqrt(square_sum);
       	 		for(int i=0;i<n;i++) B[i][j]=(B[i][j]/square_sum); // Eq. (19).
        	}

        	// Calculate the weighted normalized decision matrix based on Eq. (20)
        	for(int j=0;j<c;j++){
            	for(int i=0;i<n;i++) B[i][j]=B[i][j]*weightMatrixFN[l][j];  // Eq. (20)
        	}


        	double[] positive=new double[c];
        	double[] negative=new double[c];
        	// Evaluate the positive and negative ideal solutions for each criterion k as 
        	// minimum and maximum value of kth column weighted normalized matrix.
        	for(int j=0;j<c;j++){
            	negative[j]=positive[j]=B[0][j];
            	for(int i=0;i<n;i++){
                	negative[j]=max(negative[j],B[i][j]);
                	positive[j]=min(positive[j],B[i][j]);
            	}
        	}

        	double[] d_pos=new double[n];
        	double[] d_neg=new double[n];
        	double[] p=new double[n];
 			// Determine the distance of each alternative from the positive and negative ideal 
 			// solutions based on Eq.(21) and (22).
        	for(int i=0;i<n;i++){
            	for(int j=0;j<c;j++){
                	d_pos[i]+=(B[i][j]-positive[j])*(B[i][j]-positive[j]); 
                	d_neg[i]+=(B[i][j]-negative[j])*(B[i][j]-negative[j]);   
            	}
            	d_pos[i]=Math.sqrt(d_pos[i]);  // Eq.(21)
            	d_neg[i]=Math.sqrt(d_neg[i]);  // Eq.(22)
                                 
        	}

        	// Compute the performance score for each alternative following Eq. (23),
        	for(int i=0;i<n;i++){
        	 	p[i]=d_neg[i]/(d_neg[i]+d_pos[i]);     // Eq.(23)    
        	}

        	// ranked in decreasing order of their performance scores 
        	//Vector<Pair <float, Integer>> vv=new Vector<Pair <float, Integer>>(n);
        	PAIR[] vv=new PAIR[n];
        	for(int i=0;i<n;i++){
        		vv[i]=new PAIR(p[i],i);
        	}
        	//sort_reverse(vv,n);
        	for (int i = 0; i < n; i++) {
                for (int j = i + 1; j < n; j++) {
                    PAIR temp = new PAIR(0.0,0);
                    if (vv[j].val < vv[i].val) {
                        // Swapping
                        temp = vv[i];
                        vv[i] = vv[j];
                        vv[j] = temp;
                    }
                }
            }
            int [] p1 = new int[n];
            PAIR temp = new PAIR(0.0,0);
        	int i1=0,j1=n-1;
        	while(i1<=j1){
        	    temp = vv[i1];
                vv[i1] = vv[j1];
                vv[j1] = temp;
                i1++;
                j1--;
        	}
        	for(int i=0;i<n;i++)
        		p1[i]=vv[i].no;
        	performanceFN[l]=p1;
    	}
    
    }
    
    public static int find(Task t, int[] pl)  //find the location of task t in f's PL 
    {
        int idx =-1;
        int n= noOfTask;
        for(int i=0;i<n;i++)
        {
            if( T[pl[i]] == t)
            {
                //idx = T[pl[i]].id;
                idx = i;
                break;
            }
        }
        return idx;
    }
    
    
    public static void sort(ArrayList<pair1> list) {
 
        list.sort((o1, o2)
                  -> o2.getid().compareTo(
                      o1.getid()));
    }
    
    
    
    public static void meto_matching(){
        int n = T.length;
        int m = noOfFN;
        int k=0,idx;
        
        Fog f;
        ArrayList<ArrayList<Integer>> pltasks = new  ArrayList<ArrayList<Integer>>(n);
        ArrayList<ArrayList<Integer>> plfog = new  ArrayList<ArrayList<Integer>>(m);
        
        ArrayList<Integer> plt = new ArrayList<Integer>(m);
        ArrayList<Integer> plf = new ArrayList<Integer>(n);
        
        
        //generating pltasks
        for(int i=0;i<n;i++)
        {
            pltasks.add(new ArrayList<Integer>(m));
            for(int j=0;j<m;j++)
            {
                //pltasks.add(new ArrayList<Integer>(m));
                //plt.add(j,performanceTask[i][j]);
                pltasks.get(i).add(j,performanceTaskmeto[i][j]);
                 //System.out.println("Hi" + plt.get(j));
            }
            //pltasks.add(i,plt);
            //plt.clear();
            //System.out.println(pltasks.get(i).get(0));
            //System.out.println("yes");
        }
        
        
        //generating plfog
        for(int i=0;i<m;i++)
        {
            plfog.add(new ArrayList<Integer>(n));
            for(int j=0;j<n;j++)
            {
                //plfog.add(new ArrayList<Integer>(n));
                //plf.add(performanceFN[i][j]);
                plfog.get(i).add(performanceFN[i][j]);
            }
            //plfog.add(i,plf);
        }
        
        //Printing pltasks
        for(int i=0;i<n;i++)
        {
            System.out.println("PLTask" +i +": ");
            for(int j=0;j<m;j++)
            {
                System.out.print(pltasks.get(i).get(j)+ ", ");
            }
            System.out.println("");
        }
        
        
        //printing plfog
        for(int i=0;i<m;i++)
        {
            System.out.println("PLFog"+i+": ");
            for(int j=0;j<n;j++)
            {
                System.out.print(+plfog.get(i).get(j)+ ", ");
            }
            System.out.println("");
        }
        
         //System.out.println("Hi");
        //System.out.println(pltasks.get(0).get(0));

        //allocating space for arraylist in fm since randomly any fog node will be acccessed
        for(int i=0;i<m;i++)
         {
                   fmmeto.add(new ArrayList<pair1>(n));
         }
         
         
         //allocating space for fog node in tm since randomly any task will be acccessed
         Fog f1 = new Fog(-1,0,0,0,0,0);
         for(int i=0;i<n;i++)
         {
             tmmeto.add(i,f1);
         }
         System.out.println(tmmeto.size());
        
        //fog node allocation process
        int alloccount =0,ifcount=0,elsecount=0,rejectcount=0;
        while(Tasksmeto.isEmpty() == false) //tasklist is not empty
        {
             //System.out.println("Hi");
            k++; //pop a task;
            //get a task from tasklist and remove it from tasklist
            Task t = Tasksmeto.get(0);
            Tasksmeto.remove(0);
           // f = F[performanceTask[k][0]];
           //System.out.println("Hi3");
            
            // get the id of the most preferrred  fog node from pltasks
           int x = pltasks.get(t.id).get(0);
            //System.out.println(x);
            
           //get the  most preferred fog node and store it in f
           f = F[pltasks.get(t.id).get(0)]; 
           //System.out.println("mpf for task " +t.id+ "is : " + f.id );
           //System.out.println("Hi4");
            
           //remove the most preferred fog node from t's PL
           pltasks.get(t.id).remove(0);
           //System.out.println("Hi5");
            
           /* for(int i=0;i<n;i++)
            {
                 pfn.add(performanceFN[f.id][i]);
            }*/
             //System.out.println("Hi1");
             
            
            //find the position of task t in f's PL
            idx = find(T[t.id] , performanceFN[f.id]);
            //System.out.println(idx);
            
            //store task t and its position as a pair p1
            pair1 p1;
            p1 = new pair1(T[t.id],idx);
             //System.out.println("Hi2");
            
            if(f.vrum > 0)       //if f has vru for t
            {
                //System.out.println("f.vru>0");
                alloccount++;
                ifcount++;
               f.vrum--;
               
               //fm.add(new ArrayList<pair1>(n));
               // At 0th row, modifing the default value to 13
               //System.out.println("al1");
               
               //allocate pair p1 to f in fm i.e allocate task t to f
               //count represents the nummber of tasks allocated to f. Intitially count is zero for each FN
               fmmeto.get(f.id).add(f.countm,p1); 
               f.countm++;
               //System.out.println("al");
               
               //sort the tasks allocated to f in fm based on preference. lower idx means higher preference
               sort(fmmeto.get(f.id));
               
               //add f to t in tm
              System.out.println("t.id =" +t.id + "tm.size() =" +tmmeto.size());
               tmmeto.set(t.id,f);
              // tm.remove(t.id);
              // tm.add(t.id,f);
                //System.out.println("Hi");
            }
            else
            {
                //System.out.println("f.vru=0");
                //System.out.println("priority indx of t*" + fm.get(f.id).get(0).id);
                
                //if there exists a task with lower preference than t
                
                if(idx < fmmeto.get(f.id).get(0).id)
                {
                    alloccount++;
                    elsecount++;
                    //System.out.println("lower priority task exists");
                    pair1 p2 = fmmeto.get(f.id).get(0);  //Removing lowest priority task t*
                    Tasksmeto.add(p2.t); //add the task to be removed (p2) to Tasks list
                    fmmeto.get(f.id).remove(0); // Remove task p2 from f in fm
                    tmmeto.remove(p2.t.id); // Remove fm from t* in tm
                    tmmeto.add(p2.t.id,f1);
                    //System.out.println("hello");
                    fmmeto.get(f.id).add(f.countm - 1,p1); //assign task t to fn f, count remains same as before, we use count -1 bcoz index starts from 0
                    sort(fmmeto.get(f.id)); //sort the tasks mappped to f according to f's preferences
                    System.out.println("t.id =" +t.id + " tm.size() = " +tmmeto.size());
                    tmmeto.set(t.id,f); // aad f to t's mapping list
                    //tm.remove(t.id);
                    //tm.add(t.id,f);
                    
                    
                }
                else
                {
                    //f rejects t...it is added again to T
                    //alloccount++;
                    rejectcount++;
                    Tasksmeto.add(t);
                    
                    
                }
            }
            
            
            
        }
        
        int j;
        //Task t = new Task(-1, 0,0,0,0,0);
        
        
        //System.out.println(fm.get(0).get(0).id);
        //System.out.println(fm.get(2).get(0).id);
        
        System.out.println("---------Final Matching--METO---------");
        for(int i=0;i<m;i++)
        {
            System.out.println(" ");
            System.out.print("Fog" + i + ": " );
            j=0;
            
            while(fmmeto.get(i).isEmpty() == false && j<fmmeto.get(i).size())
            {
                System.out.print(", " + fmmeto.get(i).get(j).t.id);
                j++;
            }
            System.out.println(".");
        }
        
        System.out.println("");
        System.out.println("---------Total Values Incurred----------");
        
        //if(k==T.length && alloccount == T.length)
        //System.out.println("ok");
        //else
        //System.out.println("k=" +k + " alloccount=" +alloccount);
        
        //System.out.println("ifcount =" +ifcount + "elsecount=" +elsecount + "rejectcount =" + rejectcount);
        
        
        double l1=0 , eiot =0, efog=0, tote=0;
        for(int i=0;i<n;i++)
        {
            Fog f2 = tmmeto.get(i);
            //System.out.println("f2.id  is " + f2.id);
            l1 = l1 + decisionMatrixTaskmeto[i][f2.id][0];
            //System.out.println("f2.id  is " + f2.id);
            //System.out.println("decisionMatrixTask[i][f2.id][0]=  "+i+f2.id +"=" + decisionMatrixTask[i][f2.id][0]+ ", l1="+l1);
            eiot = eiot + (decisionMatrixTaskmeto[i][f2.id][1])/1000;
            efog = efog + decisionMatrixFN[f2.id][i][0];
            tote = tote + eiot + efog;
            
        }
        
        
        
        System.out.println("total latency:-"+l1);
        System.out.println("Average Latency per task:- " + l1/n + "sec");
        System.out.println("Total System Energy Consumption:-" + tote + "J" );
        System.out.println("Outages");
        
        meto_tote = tote;
        meto_avgl = l1/n;
        
        int o =0; //outage
        double time;
        for(int i=0;i<n;i++)
        {
            Fog f3 = tmmeto.get(i);
            time = decisionMatrixTaskmeto[i][f3.id][0];
            if(T[i].d < time)
            {
                System.out.println(i+", ");
                o++;//outages
            }
        }
        meto_outages = o;
        System.out.println("Total outages:-" +o);
        
        
    }
    
    
    
    
    
    
    public static void matching1(){
        int n = T.length;
        int m = noOfFN;
        int k=0,idx;
        
        Fog f;
        ArrayList<ArrayList<Integer>> pltasks = new  ArrayList<ArrayList<Integer>>(n);
        ArrayList<ArrayList<Integer>> plfog = new  ArrayList<ArrayList<Integer>>(m);
        
        ArrayList<Integer> plt = new ArrayList<Integer>(m);
        ArrayList<Integer> plf = new ArrayList<Integer>(n);
        
        
        //generating pltasks
        for(int i=0;i<n;i++)
        {
            pltasks.add(new ArrayList<Integer>(m));
            for(int j=0;j<m;j++)
            {
                //pltasks.add(new ArrayList<Integer>(m));
                //plt.add(j,performanceTask[i][j]);
                pltasks.get(i).add(j,performanceTask[i][j]);
                 //System.out.println("Hi" + plt.get(j));
            }
            //pltasks.add(i,plt);
            //plt.clear();
            //System.out.println(pltasks.get(i).get(0));
            //System.out.println("yes");
        }
        
        
        //generating plfog
        for(int i=0;i<m;i++)
        {
            plfog.add(new ArrayList<Integer>(n));
            for(int j=0;j<n;j++)
            {
                //plfog.add(new ArrayList<Integer>(n));
                //plf.add(performanceFN[i][j]);
                plfog.get(i).add(performanceFN[i][j]);
            }
            //plfog.add(i,plf);
        }
        
        //Printing pltasks
        for(int i=0;i<n;i++)
        {
            System.out.println("PLTask" +i +": ");
            for(int j=0;j<m;j++)
            {
                System.out.print(pltasks.get(i).get(j)+ ", ");
            }
            System.out.println("");
        }
        
        
        //printing plfog
        for(int i=0;i<m;i++)
        {
            System.out.println("PLFog"+i+": ");
            for(int j=0;j<n;j++)
            {
                System.out.print(+plfog.get(i).get(j)+ ", ");
            }
            System.out.println("");
        }
        
         //System.out.println("Hi");
        //System.out.println(pltasks.get(0).get(0));

        //allocating space for arraylist in fm since randomly any fog node will be acccessed
        for(int i=0;i<m;i++)
         {
                   fm.add(new ArrayList<pair1>(n));
         }
         
         
         //allocating space for fog node in tm since randomly any task will be acccessed
         Fog f1 = new Fog(-1,0,0,0,0,0);
         for(int i=0;i<n;i++)
         {
             tm.add(i,f1);
         }
         System.out.println(tm.size());
        
        //fog node allocation process
        int alloccount =0,ifcount=0,elsecount=0,rejectcount=0;
        while(Tasks.isEmpty() == false) //tasklist is not empty
        {
             //System.out.println("Hi");
            k++; //pop a task;
            //get a task from tasklist and remove it from tasklist
            Task t = Tasks.get(0);
            Tasks.remove(0);
           // f = F[performanceTask[k][0]];
           //System.out.println("Hi3");
            
            // get the id of the most preferrred  fog node from pltasks
           int x = pltasks.get(t.id).get(0);
            //System.out.println(x);
            
           //get the  most preferred fog node and store it in f
           f = F[pltasks.get(t.id).get(0)]; 
           //System.out.println("mpf for task " +t.id+ "is : " + f.id );
           //System.out.println("Hi4");
            
           //remove the most preferred fog node from t's PL
           pltasks.get(t.id).remove(0);
           //System.out.println("Hi5");
            
           /* for(int i=0;i<n;i++)
            {
                 pfn.add(performanceFN[f.id][i]);
            }*/
             //System.out.println("Hi1");
             
            
            //find the position of task t in f's PL
            idx = find(T[t.id] , performanceFN[f.id]);
            //System.out.println(idx);
            
            //store task t and its position as a pair p1
            pair1 p1;
            p1 = new pair1(T[t.id],idx);
             //System.out.println("Hi2");
            
            if(f.vru > 0)       //if f has vru for t
            {
                //System.out.println("f.vru>0");
                alloccount++;
                ifcount++;
               f.vru--;
               
               //fm.add(new ArrayList<pair1>(n));
               // At 0th row, modifing the default value to 13
               //System.out.println("al1");
               
               //allocate pair p1 to f in fm i.e allocate task t to f
               //count represents the nummber of tasks allocated to f. Intitially count is zero for each FN
               fm.get(f.id).add(f.count,p1); 
               f.count++;
               //System.out.println("al");
               
               //sort the tasks allocated to f in fm based on preference. lower idx means higher preference
               sort(fm.get(f.id));
               
               //add f to t in tm
              System.out.println("t.id =" +t.id + "tm.size() =" +tm.size());
               tm.set(t.id,f);
              // tm.remove(t.id);
              // tm.add(t.id,f);
                //System.out.println("Hi");
            }
            else
            {
                //System.out.println("f.vru=0");
                //System.out.println("priority indx of t*" + fm.get(f.id).get(0).id);
                
                //if there exists a task with lower preference than t
                
                if(idx < fm.get(f.id).get(0).id)
                {
                    alloccount++;
                    elsecount++;
                    //System.out.println("lower priority task exists");
                    pair1 p2 = fm.get(f.id).get(0);  //Removing lowest priority task t*
                    Tasks.add(p2.t); //add the task to be removed (p2) to Tasks list
                    fm.get(f.id).remove(0); // Remove task p2 from f in fm
                    tm.remove(p2.t.id); // Remove fm from t* in tm
                    tm.add(p2.t.id,f1);
                    //System.out.println("hello");
                    fm.get(f.id).add(f.count - 1,p1); //assign task t to fn f, count remains same as before, we use count -1 bcoz index starts from 0
                    sort(fm.get(f.id)); //sort the tasks mappped to f according to f's preferences
                    System.out.println("t.id =" +t.id + " tm.size() = " +tm.size());
                    tm.set(t.id,f); // aad f to t's mapping list
                    //tm.remove(t.id);
                    //tm.add(t.id,f);
                    
                    
                }
                else
                {
                    //f rejects t...it is added again to T
                    //alloccount++;
                    rejectcount++;
                    Tasks.add(t);
                    
                    
                }
            }
            
            
            
        }
        
        int j;
        //Task t = new Task(-1, 0,0,0,0,0);
        
        
        //System.out.println(fm.get(0).get(0).id);
        //System.out.println(fm.get(2).get(0).id);
        
        System.out.println("---------Final Matching--METO---------");
        for(int i=0;i<m;i++)
        {
            System.out.println(" ");
            System.out.print("Fog" + i + ": " );
            j=0;
            
            while(fm.get(i).isEmpty() == false && j<fm.get(i).size())
            {
                System.out.print(", " + fm.get(i).get(j).t.id);
                j++;
            }
            System.out.println(".");
        }
        
        System.out.println("");
        System.out.println("---------Total Values Incurred----------");
        
        //if(k==T.length && alloccount == T.length)
        //System.out.println("ok");
        //else
        //System.out.println("k=" +k + " alloccount=" +alloccount);
        
        //System.out.println("ifcount =" +ifcount + "elsecount=" +elsecount + "rejectcount =" + rejectcount);
        
        
        double l1=0 , eiot =0, efog=0, tote=0;
        for(int i=0;i<n;i++)
        {
            Fog f2 = tm.get(i);
            //System.out.println("f2.id  is " + f2.id);
            l1 = l1 + decisionMatrixTask[i][f2.id][0];
            //System.out.println("f2.id  is " + f2.id);
            //System.out.println("decisionMatrixTask[i][f2.id][0]=  "+i+f2.id +"=" + decisionMatrixTask[i][f2.id][0]+ ", l1="+l1);
            eiot = eiot + (decisionMatrixTask[i][f2.id][1])/1000;
            efog = efog + decisionMatrixFN[f2.id][i][0];
            tote = tote + eiot + efog;
            
        }
        
        System.out.println("total latency:-"+l1);
        System.out.println("Average Latency per task:- " + l1/n + "sec");
        System.out.println("Total System Energy Consumption:-" + tote + "J" );
        System.out.println("Outages");
        
        mapiot_tote = tote;
        mapiot_avgl = l1/n;
        
        int o =0; //outage
        double time;
        for(int i=0;i<n;i++)
        {
            Fog f3 = tm.get(i);
            time = decisionMatrixTask[i][f3.id][0];
            if(T[i].d < time)
            {
                System.out.println(i+", ");
                o++;//outages
            }
        }
        mapiot_outage = o;
        System.out.println("Total outages:-" +o);
        
        
    }



    
    public static void vmap_matching()
    {
        
        int n = T.length;
        int m = noOfFN;
        int idx;
        double dealloc =0, temp=0;
        Fog f;
        ArrayList<ArrayList<Integer>> pltasks = new  ArrayList<ArrayList<Integer>>(n);
        ArrayList<ArrayList<Integer>> plfog = new  ArrayList<ArrayList<Integer>>(m);
        
        ArrayList<Integer> plt = new ArrayList<Integer>(m);
        ArrayList<Integer> plf = new ArrayList<Integer>(n);
        
        
        //generating pltasks
        for(int i=0;i<n;i++)
        {
            pltasks.add(new ArrayList<Integer>(m));
            for(int j=0;j<m;j++)
            {
                //pltasks.add(new ArrayList<Integer>(m));
                //plt.add(j,performanceTask[i][j]);
                pltasks.get(i).add(j,performanceTaskvmap[i][j]);
                 //System.out.println("Hi" + plt.get(j));
            }
            //pltasks.add(i,plt);
            //plt.clear();
            //System.out.println(pltasks.get(i).get(0));
            //System.out.println("yes");
        }
        
        
        //generating plfog
        for(int i=0;i<m;i++)
        {
            plfog.add(new ArrayList<Integer>(n));
            for(int j=0;j<n;j++)
            {
                //plfog.add(new ArrayList<Integer>(n));
                //plf.add(performanceFN[i][j]);
                plfog.get(i).add(performanceFNvmap[i][j]);
            }
            //plfog.add(i,plf);
        }
        
        //Printing pltasks
        for(int i=0;i<n;i++)
        {
            System.out.println("PLTask" +i +": ");
            for(int j=0;j<m;j++)
            {
                System.out.print(pltasks.get(i).get(j)+ ", ");
            }
            System.out.println("");
        }
        
        
        //printing plfog
        for(int i=0;i<m;i++)
        {
            System.out.println("PLFog"+i+": ");
            for(int j=0;j<n;j++)
            {
                System.out.print(+plfog.get(i).get(j)+ ", ");
            }
            System.out.println("");
        }
        
         //System.out.println("Hi");
        //System.out.println(pltasks.get(0).get(0));

        //allocating space for arraylist in fm since randomly any fog node will be acccessed
        for(int i=0;i<m;i++)
         {
                   fmvmap.add(new ArrayList<pair1>(n));
         }
         
         
         //allocating space for fog node in tm since randomly any task will be acccessed
         Fog f1 = new Fog(-1,0,0,0,0,0);
         for(int i=0;i<n;i++)
         {
             tmvmap.add(i,f1);
         }
         System.out.println(tmvmap.size());
         System.out.println("Hi");
        
        
        
        
        int y=0,k=0;
        Task t = new Task(-1,0,0,0,0,0);
        Task t1 = new Task(-2,0,0,0,0,0);
        f = new Fog(-1,-1,-1,-1,-1,-1);
        
        System.out.println("Hi");
        
       // Task t1 = new Task(-1,0,0,0,0,0);
        while(Tasksvmap.isEmpty() == false)
        {
            y=0;
            //Task t = new Task(-1,-1,-1,-1,-1,-1);
            t = Tasksvmap.get(0);
            Tasksvmap.remove(0);
            //Fog f = new Fog(-1,-1,-1,-1,-1);
            f = F[pltasks.get(t.id).get(0)];
            pltasks.get(t.id).remove(0);
            idx = find(t, performanceFNvmap[f.id]);
            pair1 p = new pair1(t,idx);
            
            
            System.out.println("Hi");
            
            
            if(compr[t.id][f.id] < f.avail_cpu)
            {
                y=1;
                fmvmap.get(f.id).add(p);
                sort(fmvmap.get(f.id));
                tmvmap.add(f);   //
                f.avail_cpu = f.avail_cpu - compr[t.id][f.id];
                
            }
            else
            {
                
                if(compr[t.id][f.id] > f.avail_cpu)
                {
                    //y=0;
                    dealloc =0;
                    temp=f.avail_cpu;
                    //find task having less priority than than
                    k=0;
                    while( k < fmvmap.get(f.id).size() )
                    {
                        if(p.id < fmvmap.get(f.id).get(k).id)
                        {
                            //y=1;
                            //k++;
                            t1 = fmvmap.get(f.id).get(k).t;
                            temp = temp + compr[t1.id][f.id];
                            dealloc = dealloc + compr[t1.id][f.id];
                            k++;
                            if(compr[t.id][f.id] < temp)
                            {
                                y=1;
                                //allocate task t to f
                                fmvmap.get(f.id).add(p);
                                sort(fm.get(f.id));
                                tmvmap.set(t.id,f);
                                f.avail_cpu = f.avail_cpu - dealloc + compr[t.id][f.id];
                                break;
                                
                            }
                        }
                    }
                    
                    if(y==1)
                    {
                        for(int i=0;i<k;i++)
                        {
                            fmvmap.get(f.id).remove(i);
                            
                        }
                    }
                }
                
            }
            
            if(y==0)
            {
                Tasksvmap.add(t);
                
            }
            
        }
    }
    
    
    

    public static void matching(){


    // Q -- Quota of each FN
   	Integer[] Q=new Integer[noOfFN];
   	for(int i=0;i<noOfFN;i++)
   	    Q[i]=200; // store some value

    // Assign  -- Assigned tasks in FN
    Integer[][] Assign=new Integer[noOfFN][noOfTask];
    	for(int i=0;i<noOfTask;i++){
    	    int[] tj=new int[noOfFN];
    	    tj=performanceTask[i];
        	int n=noOfFN;
        	for(int x=0;x<n;x++){
        		// fiÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã¢â‚¬Â¹ÃƒÂ¢Ã¢â€šÂ¬Ã‚Â ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬ÃƒÂ¢Ã¢â€šÂ¬Ã‚Â= highest ranked FN in P(tj) to which tj has not proposed yet
        		// Send proposal to fiÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã¢â‚¬Â¹ÃƒÂ¢Ã¢â€šÂ¬Ã‚Â ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬ÃƒÂ¢Ã¢â€šÂ¬Ã‚Â.
            	int fi=(int)tj[x];
            	if(Q[fi]>0){    // if QiÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã¢â‚¬Â¹ÃƒÂ¢Ã¢â€šÂ¬Ã‚Â ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬ÃƒÂ¢Ã¢â€šÂ¬Ã‚Â > 0 then
                	Assign[fi][x]=i+100;
                	System.out.println("Task "+i+" is being assigned to Fog Device "+fi);
                	Q[fi]=Q[fi]-1;
                	break;
            	}
            	else{
                	// Reject the assignment request;
           	}
        }
    	
    	    
    	}
    	


	}
	
	
	
	private static void print1(){
	    System.out.println("Task decisionMatrix");
	    for(int i=0;i<noOfTask;i++){
	        System.out.println("Task"+i+" -----------");
	        for(int j=0;j<noOfFN;j++){
	            System.out.println(decisionMatrixTask[i][j][0]+" "+decisionMatrixTask[i][j][1]);
	        }
	        
	    }
	    
	    System.out.println("FN decisionMatrix");
	    for(int i=0;i<noOfFN;i++){
	        System.out.println("FN"+i+" -----------");
	        for(int j=0;j<noOfTask;j++){
	            System.out.println(decisionMatrixFN[i][j][0]+" "+decisionMatrixFN[i][j][1]);
	        }
	        
	    }
	}
	
	private static void print2(){
	    System.out.println("Task weightMatrixTask");
	    for(int i=0;i<noOfTask;i++){
	        System.out.println("Task"+i+" -----------");
	        System.out.println(weightMatrixTask[i][0]+" "+weightMatrixTask[i][1]);
	        
	        
	    }
	    
	    System.out.println("FN weightMatrixFN");
	    for(int i=0;i<noOfFN;i++){
	        System.out.println("FN"+i+" -----------");
	       
	        System.out.println(weightMatrixFN[i][0]+" "+weightMatrixFN[i][1]);
	        
	        
	    }
	}
	
	
	private static void print3(){
	    System.out.println("Task performanceTask");
	    for(int i=0;i<noOfTask;i++){
	        System.out.println("Task"+i+" -----------");
	        for(int j=0;j<noOfFN;j++)
	            System.out.println(performanceTask[i][j]);
	      
	        
	    }
	    
	    System.out.println("FN performanceFN");
	    for(int i=0;i<noOfFN;i++){
	        System.out.println("FN"+i+" -----------");
	       
	        for(int j=0;j<noOfTask;j++)
	            System.out.println(performanceFN[i][j]);
	        
	    }
	}



	
	
	
	
	
	
	
	
	private static double generateRandom(double min, double max ) {

    // find diff
    double difference = max - min;

    // generate random number 
    double rand = Math.random();

    // multiply with difference 
    rand = Math.floor( rand * difference);

    // add with min value 
    rand = rand + min;

    return rand;

}

private static int generateRandomi(int min, int max ) {

    // find diff
    int difference = max - min;

    // generate random number 
    int rand = (int)Math.round(Math.random());

    // multiply with difference 
    rand = (int)Math.round(Math.floor( rand * difference));

    // add with min value 
    rand = rand + min;

    return rand;
	
}


}