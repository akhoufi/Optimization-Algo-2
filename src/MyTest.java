import java.util.StringTokenizer;

public class MyTest {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		String property = System.getProperty("java.library.path");
		StringTokenizer parser = new StringTokenizer(property, ";");
		while (parser.hasMoreTokens()) {
		    System.err.println(parser.nextToken());
		    }
	}

}
