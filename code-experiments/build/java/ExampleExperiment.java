import java.util.Random;
import java.util.TreeMap;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
/**
 * An example of benchmarking random search on a COCO suite. 
 *
 * Set the parameter BUDGET_MULTIPLIER to suit your needs.
 */
public class ExampleExperiment {

	/**
	 * The maximal budget for evaluations done by an optimization algorithm equals 
	 * dimension * BUDGET_MULTIPLIER.
	 * Increase the budget multiplier value gradually to see how it affects the runtime.
	 */
	public static int BUDGET_MULTIPLIER = 2;
	
	/**
	 * The maximal number of independent restarts allowed for an algorithm that restarts itself. 
	 */
	public static final int INDEPENDENT_RESTARTS = 10000;
	
	/**
	 * The random seed. Change if needed.
	 */
	public static final long RANDOM_SEED = 0xdeadbeef;

	/**
	 * The problem to be optimized (needed in order to simplify the interface between the optimization
	 * algorithm and the COCO platform).
	 */
	public static Problem PROBLEM;
	
	/**
	 * Interface for function evaluation.
	 */
	public interface Function {
		double[] evaluate(double[] x);
    }

	/**
	 * Evaluate the static PROBLEM.
	 */
    public static final Function evaluateFunction = new Function() {
    	public double[] evaluate(double[] x) {
    		return PROBLEM.evaluateFunction(x);
    	}
    };

	/**
	 * The main method initializes the random number generator and calls the example experiment on the
	 * bi-objective suite.
	 */
	public static void main(String[] args) {
		
		Random randomGenerator = new Random(RANDOM_SEED);
		
		if(args.length > 0){
			BUDGET_MULTIPLIER = Integer.valueOf(args[0]);
		}
		
		/* Change the log level to "warning" to get less output */
		CocoJNI.cocoSetLogLevel("info");

		System.out.println("Running the example experiment... (might take time, be patient)");
		System.out.flush();
		
		/* Call the example experiment 
		exampleExperiment("bbob-biobj", "bbob-biobj", randomGenerator);*/

		/* Uncomment the line below to run the same example experiment on the bbob suite */
	  	exampleExperiment("bbob", "bbob", randomGenerator); 

		System.out.println("Done!");
		System.out.flush();

		return;
	}
	
	/**
	 * A simple example of benchmarking random search on a suite with instances from 2016.
	 *
	 * @param suiteName Name of the suite (use "bbob" for the single-objective and "bbob-biobj" for the
	 * bi-objective suite).
	 * @param observerName Name of the observer (use "bbob" for the single-objective and "bbob-biobj" for the
	 * bi-objective observer).
	 * @param randomGenerator The random number generator.
	 */
	public static void exampleExperiment(String suiteName, String observerName, Random randomGenerator) {
		try {

			/* Set some options for the observer. See documentation for other options. */
			final String observerOptions = 
					  "result_folder: RS_on_" + suiteName + " " 
					+ "algorithm_name: RS "
					+ "algorithm_info: \"A simple random search algorithm\"";

			/* Initialize the suite and observer */
			Suite suite = new Suite(suiteName, "year: 2016", "dimensions: 2,3,5,10,20,40");
			Observer observer = new Observer(observerName, observerOptions, (long)BUDGET_MULTIPLIER);
			Benchmark benchmark = new Benchmark(suite, observer);

			/* Initialize timing */
			Timing timing = new Timing();
			
			/* Iterate over all problems in the suite */
			while ((PROBLEM = benchmark.getNextProblem()) != null) {

				int dimension = PROBLEM.getDimension();

				/* Run the algorithm at least once */
				for (int run = 1; run <= 1 + INDEPENDENT_RESTARTS; run++) {

					long evaluationsDone = PROBLEM.getEvaluations();
					long evaluationsRemaining = (long) (dimension * BUDGET_MULTIPLIER) - evaluationsDone;

					/* Break the loop if the target was hit or there are no more remaining evaluations */
					if (PROBLEM.isFinalTargetHit() || (evaluationsRemaining <= 0))
						break;

					/* Call the optimization algorithm for the remaining number of evaluations */
					myRandomSearch(evaluateFunction,
							       dimension,
							       PROBLEM.getNumberOfObjectives(),
							       PROBLEM.getSmallestValuesOfInterest(),
							       PROBLEM.getLargestValuesOfInterest(),
							       evaluationsRemaining,
							       randomGenerator);

					/* Break the loop if the algorithm performed no evaluations or an unexpected thing happened */
					if (PROBLEM.getEvaluations() == evaluationsDone) {
						System.out.println("WARNING: Budget has not been exhausted (" + evaluationsDone + "/"
								+ dimension * BUDGET_MULTIPLIER + " evaluations done)!\n");
						break;
					} else if (PROBLEM.getEvaluations() < evaluationsDone)
						System.out.println("ERROR: Something unexpected happened - function evaluations were decreased!");
				}

				/* Keep track of time */
				timing.timeProblem(PROBLEM);
			}

			/* Output the timing data */
			timing.output();

			benchmark.finalizeBenchmark();

		} catch (Exception e) {
			System.err.println(e.toString());
		}
	}

	/** 
	 * A simple random search algorithm that can be used for single- as well as multi-objective 
	 * optimization.
	 */
	public static void myRandomSearch(Function f, 
			                          int dimension, 
			                          int numberOfObjectives, 
			                          double[] lowerBounds,
			                          double[] upperBounds, 
			                          long maxBudget, 
			                          Random randomGenerator) {

		double[] x = new double[dimension];
		double[] sigma = new double[dimension];
		double[] y = new double[numberOfObjectives];
		double range;
		
		// Declaration of parameter
		// Keeping the lambda low ==> keeping the no. of evaluations low ==> more helpful for search I guess
		int lambda = 5*dimension; long mu = Math.floorDiv(lambda, 4) ; double tau = 1/Math.sqrt(dimension); double tau_i = Math.sqrt(tau);
		
		

	    /* Construct x as a random point between the lower and upper bounds. Also constructing initial sigma */
		for (int j = 0; j < dimension; j++) {
			range = upperBounds[j] - lowerBounds[j];
			x[j] = lowerBounds[j] + randomGenerator.nextDouble() * range;
			sigma[j] = (upperBounds[j] - lowerBounds[j])/20.0;	// 1/4th of the spread in that dimension
		}
		
		// Budget is defined as the no. of f-evals, since each time we create new population, we do lambda evals, so modifying myBudget accordingly 
		long myBudget = Math.floorDiv(maxBudget, lambda);
		if (myBudget < 1){
			lambda = (int) maxBudget;		//No. of new points generated = maxBudget now, which is now smaller than original lambda
			myBudget = 1;
		}

		//This Loop is the replacement of line-2 ("While not happy")		
		for (int i = 0; i < myBudget; i++) {
			
			
			double[] epsilon = new double[lambda];
			double[][] epsilon_vec = new double[lambda][dimension];
			double[][] points_vec = new double[lambda][dimension];
			double[][] change_vec = new double[lambda][dimension];
			double[][] sigma_vec = new double[lambda][dimension];

			// Inner-loop for generating lambda new-elements starts
			for (int k=0; k < lambda; k++){
				epsilon[k] = tau* randomGenerator.nextGaussian();
				
				// Generating epsilon_k
				for (int j=0; j < dimension; j++){
					epsilon_vec[k][j] = tau_i*randomGenerator.nextGaussian();
				}
				
				// Generating change_k
				for (int j=0; j < dimension; j++){
					change_vec[k][j] = randomGenerator.nextGaussian();
				}

				// Generating sigma_k
				for (int j=0; j < dimension; j++){
					sigma_vec[k][j] = sigma[j]*Math.exp(epsilon_vec[k][j])*Math.exp(epsilon[k]);
				}
				
				// Generating x_k, note that no bounds check is performed here for the new points
				for (int j=0; j < dimension; j++){
					points_vec[k][j] = x[j]+sigma_vec[k][j]*change_vec[k][j];
				}
				
			}
			// Inner-loop for generating lambda new-elements ends
			
			
			// Calculating fitness values at new points
			double[] f_vals = new double[lambda];
			for (int k=0; k < lambda; k++){
				
				double[] point = new double[dimension];
				for(int j=0; j < dimension; j++){
					point[j] = points_vec[k][j];
				}
			
				y = f.evaluate(point);
				f_vals[k] = y[0];	//Adding zero-index as our algorithm is designed for single-objective problems
			}
			
			//Sorting ascending order
			List<Integer> indices = dirtySort(f_vals);
			
			//Calculating new sigma and x
			double[] temp_x = new double[dimension];Arrays.fill(temp_x, 0.0);
			double[] temp_sigma = new double[dimension];Arrays.fill(temp_sigma, 0.0);
			

			for (int k=0; k < mu; k++){
				int index = indices.get(k);
				for (int j=0; j < dimension; j++){
					temp_x[j] += (points_vec[index][j])/mu;
					temp_sigma[j] += (sigma_vec[index][j])/mu;
				}
			}
			
			for(int j=0;j<dimension;j++){
				x[j] = temp_x[j];
				sigma[j] = temp_sigma[j];
			}
	
		}
		

	}
	
	public static List<Integer> dirtySort(double[] array){
		TreeMap<Double, List<Integer>> map = new TreeMap<Double, List<Integer>>();
		
		for(int i = 0; i < array.length; i++) {
		    List<Integer> ind = map.get(array[i]);
		    if(ind == null){
		        ind = new ArrayList<Integer>();
		        map.put(array[i], ind);
		    }
		    ind.add(i);
		}

		// Now flatten the list
		List<Integer> indices = new ArrayList<Integer>();
		for(List<Integer> arr : map.values()) {
		    indices.addAll(arr);
		}
		
		return indices;
	}
}
