import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.Stack;

class Edmond {
    static String filename;
    static ArrayList<Integer> verticesCount = new ArrayList<>();        // Store the number of vertices in each graph
    static ArrayList<ArrayList<int[]>> graphs_read = new ArrayList<>(); // Store the graphs read from the file
    /**
     * Reverses the given directed graph to obtain edges in the opposite direction (from end to start nodes).
     * @param G The input directed graph represented as a map from source nodes to their outgoing edges and weights.
     * @return The reversed graph represented as a map from destination nodes to their incoming edges and weights.
     */
    public static Map<Integer, Map<Integer, Integer>> reverseGraph(Map<Integer, Map<Integer, Integer>> G) {
        // Create a new HashMap to store the reversed graph.
        Map<Integer, Map<Integer, Integer>> g = new HashMap<>();

        // Iterate through the source nodes in the input graph.
        for (int src : G.keySet()) {
            // Iterate through the destination nodes and their corresponding weights in the current source node.
            for (int dst : G.get(src).keySet()) {
                // If the destination node is not in the reversed graph, create a new entry for it.
                if (!g.containsKey(dst)) {
                    g.put(dst, new HashMap<>());
                }
                // Reverse the edge by putting the source node as the destination and its weight.
                // This effectively reverses the direction of the edge.
                g.get(dst).put(src, G.get(src).get(dst));   
            }
        }

        return g;
    }

    /**
     * Builds the minimum graph from the given reversed graph.
     * @param rg The reversed graph.
     * @param root The root node of the graph.
     * @return The minimum graph represented as a map from destination nodes to their incoming edges and weights.
     */
    public static Map<Integer, Map<Integer, Integer>> buildMin(Map<Integer, Map<Integer, Integer>> rg, int root) {
        Map<Integer, Map<Integer, Integer>> mg = new HashMap<>();
        // Iterate through the destination nodes in the reversed graph.
        for (int dst : rg.keySet()) {
            if (dst == root) {
                continue;
            }
            int minInd = -1;
            int minValue = Integer.MAX_VALUE;
            // Iterate through the source nodes and their corresponding weights in the current destination node.
            // Find the minimum-weight edge (u, v) such that u is in S and v is not in S.
            for (int src : rg.get(dst).keySet()) {
                if (rg.get(dst).get(src) <= minValue) {
                    minInd = src;
                    minValue = rg.get(dst).get(src);
                }
            }
            // Add the minimum-weight edge (u, v) to the minimum graph.
            Map<Integer, Integer> minMap = new HashMap<>();
            minMap.put(minInd, minValue);
            mg.put(dst, minMap);
        }
        return mg;
    }
    /**
     * Find the cycle in the graph
     * @param mg The minimum graph
     * @return The cycle in the graph
     */
    public static List<Integer> findCycle(Map<Integer, Map<Integer, Integer>> mg) {
        // Iterate through the source nodes in the minimum graph.
        for (int start : mg.keySet()) {
            List<Integer> visited = new ArrayList<>();
            Stack<Integer> stack = new Stack<>();
            stack.push(start);
            // DFS to find the cycle in the graph
            // If the node is visited, then it is in the cycle and return the cycle
            // If the node is not visited, then add it to the visited list and push its neighbors to the stack
            // If the stack is empty, then there is no cycle in the graph
            while (!stack.isEmpty()) {
                int n = stack.pop();
                if (visited.contains(n)) {
                    List<Integer> C = new ArrayList<>();
                    while (!C.contains(n)) {
                        C.add(n);
                        n = new ArrayList<>(mg.get(n).keySet()).get(0);
                    }
                    return C;
                }
                visited.add(n);
                if (mg.containsKey(n)) {
                    stack.addAll(mg.get(n).keySet());
                }
            }
        }
        return null;
    }


    /**
     * This method implements Edmond's algorithm to find the minimum spanning arborescence of a directed graph.
     * It takes in a directed graph represented as a Map of Maps, where the outer map represents the source node and the inner map represents the destination node and the weight of the edge.
     * It also takes in the root node from which the minimum spanning arborescence is to be found.
     * The algorithm first reverses the graph and includes the root node in the reversed graph.
     * It then builds the minimum graph and finds the cycle in the graph.
     * If there is no cycle in the graph, then the minimum graph is the minimum arborescence.
     * If there is a cycle in the graph, then it finds the minimum weight edge in the cycle and removes it, and then runs the algorithm again.
     * The method returns the minimum spanning arborescence as a Map of Maps, where the outer map represents the source node and the inner map represents the destination node and the weight of the edge.
     *
     * @param G    the directed graph represented as a Map of Maps
     * @param root the root node from which the minimum spanning arborescence is to be found
     * @return the minimum spanning arborescence as a Map of Maps
     */
    public static Map<Integer, Map<Integer, Integer>> edmondMin(Map<Integer, Map<Integer, Integer>> G, Integer root) {
        // Reverse the graph.
        Map<Integer, Map<Integer, Integer>> rg = reverseGraph(G);
        // Including the root node in the reversed graph
        rg.put(root, new HashMap<>());

        // Build the minimum graph
        Map<Integer, Map<Integer, Integer>> mg = buildMin(rg, root);

        // Find the cycle in the graph
        List<Integer> C = findCycle(mg);

        // If there is no cycle in the graph, then te minimum graph is the minimum arborescence
        if (C == null) {
            return reverseGraph(mg);
        }

        // If there is a cycle in the graph, then find the minimum weight edge in the cycle and remove it
        // Then run the algorithm again

        // New vertex vc having key = max(G) + 1 to replace the cycle in the graph having a unique key
        int vc = Collections.max(G.keySet()) + 1;

        // Remove the cycle in the graph and add the new vertex vc
        List<Integer> all_nodes = new ArrayList<>(G.keySet());
        List<Integer> V_prime = new ArrayList<>(new HashSet<>(all_nodes));
        V_prime.removeAll(C);
        V_prime.add(vc);

        // The new graph now is G_prime
        Map<Integer, Map<Integer, Integer>> G_prime = new HashMap<>();
        Map<Integer, Integer> vc_in_idx = new HashMap<>();
        Map<Integer, Integer> vc_out_idx = new HashMap<>();

        // Iterate through the nodes in the graph
        for (int u : all_nodes) {
            // For each node, iterate through its neighbors 
            for (int v : G.get(u).keySet()) {
                // Case 1: If u is not in the cycle and v is in the cycle
                if (!C.contains(u) && C.contains(v)) {
                    // If u is not in G_prime, add it with an empty map
                    if (!G_prime.containsKey(u)) {
                        G_prime.put(u, new HashMap<>());
                    }
                    // Calculate the weight of the edge (u, vc)
                    int w = G.get(u).get(v) - new ArrayList<>(mg.get(v).values()).get(0);
                    // If (u, vc) is not in G_prime or the weight of (u, vc) is greater than w, update G_prime and vc_in_idx
                    if (!G_prime.get(u).containsKey(vc) || (G_prime.get(u).containsKey(vc) && w < G_prime.get(u).get(vc))) {
                        G_prime.get(u).put(vc, w);
                        vc_in_idx.put(u, v);
                    }
                }
                // Case 2: If u is in the cycle and v is not in the cycle
                else if (C.contains(u) && !C.contains(v)) {
                    // If vc is not in G_prime, add it with an empty map
                    if (!G_prime.containsKey(vc)) {
                        G_prime.put(vc, new HashMap<>());
                    }
                    // Calculate the weight of the edge (vc, v)
                    int w = G.get(u).get(v);
                    // If (vc, v) is not in G_prime or the weight of (vc, v) is greater than w, update G_prime and vc_out_idx
                    if (!G_prime.get(vc).containsKey(v) || (G_prime.get(vc).containsKey(v) && w < G_prime.get(vc).get(v))) {
                        G_prime.get(vc).put(v, w);
                        vc_out_idx.put(v, u);
                    }
                } 
                // Case 3: If neither u nor v is in the cycle
                else if (!C.contains(u) && !C.contains(v)) {
                    // If u is not in G_prime, add it with an empty map
                    if (!G_prime.containsKey(u)) {
                        G_prime.put(u, new HashMap<>());
                    }
                    // Add the edge (u, v) to G_prime
                    G_prime.get(u).put(v, G.get(u).get(v));
                }
            }
        }

        // Recursively run the algorithm on the new graph G_prime and get the minimum arborescence A
        Map<Integer, Map<Integer, Integer>> A = edmondMin(G_prime, root);

        // All the nodes in A 
        List<Integer> all_nodes_A = new ArrayList<>(A.keySet());

        int orig_in = -1;
        // Iterate through the nodes in A
        for (int src : all_nodes_A) {
            
            if (src == vc) {
                // Iterates over all the nodes that have an incoming edge from vc and adds them to the graph A
                for (int node_in : A.get(src).keySet()) {
                    int orig_out = vc_out_idx.get(node_in);
                    // If the node is not already in A, it creates a new HashMap for it
                    if (!A.containsKey(orig_out)) {
                        A.put(orig_out, new HashMap<>());
                    }
                    // Adds the incoming edge from vc to the node in A with the same weight as in the original graph G
                    A.get(orig_out).put(node_in, G.get(orig_out).get(node_in));
                }
            }
            // If the node is not equal to vc 
            else {
                // Creates a deep copy of the nodes in A that have an outgoing edge from the current node
                List<Integer> deep_copy = new ArrayList<>(A.get(src).keySet());
                // Iterates over each node in the deep copy
                for (int dst : deep_copy) {
                    // for the vc node
                    if (dst == vc) {
                        // add the incoming edge from the current node to vc to A wit the same weight
                        orig_in = vc_in_idx.get(src);
                        A.get(src).put(orig_in, G.get(src).get(orig_in));
                        A.get(src).remove(dst);
                    }
                }
            }
        }

        try {
            A.remove(vc);
        } catch (Exception e) {
            // do nothing
        }

        // Iterate through the nodes in the cycle
        for (int node : C) {
            if (node != orig_in) {
                int src = new ArrayList<>(mg.get(node).keySet()).get(0);
                if (!A.containsKey(src)) {
                    A.put(src, new HashMap<>());
                }
                A.get(src).put(node, mg.get(node).get(src));
            }
        }

        return A;
    }
    public static void readFile() {
        try {
            BufferedReader br = new BufferedReader(new FileReader(filename));
            String line;
            ArrayList<int[]> edges = new ArrayList<>();
            try {
                while ((line = br.readLine()) != null) {
                    if (line.contains("G")) {
                        verticesCount.add(Integer.parseInt(line.split(" ")[4]));
                    } else if (line.contains("{")) {
                        continue;
                    } else if (line.contains("(")) {
                        int[] edge = new int[3];
                        int idx = 0;
                        while (line.charAt(idx) != '(') {
                            idx++;
                        }
                        Scanner scanner = new Scanner(line.substring(idx + 1, line.length() - 1));
                        scanner.useDelimiter("[^\\d.]+");
                        edge[0] = scanner.nextInt();
                        edge[1] = scanner.nextInt();
                        edge[2] = (int) (scanner.nextDouble() * 100000);
                        scanner.close();
                        edges.add(edge);
                    } else if (line.contains("}")) {
                        break;
                    } else if (line.contains("--")) {
                        graphs_read.add(edges);
                        edges = new ArrayList<>();
                    }
                }
                br.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
            
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    public static void main(String[] args) {
        
        // From args parameter read the file
        try {
            filename = args[0];
        } catch (Exception e) {
            filename = "dwg1.txt";
        }
        // Reading the file
        readFile();

        // Minimum Arborescences in dwg1.txt
        System.out.println("Minimum Arborescences in " + filename);

        int i = 1;
        for (ArrayList<int[]> edges : graphs_read) {
            System.out.println("G" + i + ":" + "|V|=" + verticesCount.get(i-1) + " -----------------------");
            i++;
            runEdmond(edges);
        }
    }
    private static void runEdmond(ArrayList<int[]> edges) {

        Map<Integer, Map<Integer, Integer>> G = new HashMap<>();
        for (int[] edge : edges) {
            int src = edge[0];
            int dst = edge[1];
            Integer w = edge[2];
            if (!G.containsKey(src)) {
                G.put(src, new HashMap<>());
            }
            G.get(src).put(dst, w);
        }

        // Record the time in ms
        long startTime = System.currentTimeMillis();
        // Run the algorithm
        Map<Integer, Map<Integer, Integer>> A = edmondMin(G, 0);
        // Record the time in ms
        long endTime = System.currentTimeMillis();

        // Sort the graph
        A = sortTheGraph(A);

        // Output the arborescence
        System.out.println("Arborescence --");
        double totalWeight = 0;
        for (int src : A.keySet()) {
            for (int dst : A.get(src).keySet()) {
                int weight = A.get(src).get(dst);
                System.out.println(src + 1 + ": (" + src + ", " + dst + ", " + String.format("%.5f", (double) weight / 100000) + ")");
                totalWeight += (double) weight / 100000;
            }
        }
        System.out.println("Total weight: " + String.format("%.5f", totalWeight) + " (" + (endTime - startTime) + " ms)");
    }
    private static Map<Integer, Map<Integer, Integer>> sortTheGraph(Map<Integer, Map<Integer, Integer>> A) {
        /*
         *  Sort the graph to the ascending order of x (leaving vertices), and if x are the same, use y (ending vertices) as the secondary key.
         */
        List<Integer> all_nodes = new ArrayList<>(A.keySet());
        Collections.sort(all_nodes);
        Map<Integer, Map<Integer, Integer>> sorted_A = new HashMap<>();
        for (int src : all_nodes) {
            List<Integer> all_dst = new ArrayList<>(A.get(src).keySet());
            Collections.sort(all_dst);
            Map<Integer, Integer> sorted_dst = new HashMap<>();
            for (int dst : all_dst) {
                sorted_dst.put(dst, A.get(src).get(dst));
            }
            sorted_A.put(src, sorted_dst);
        }
        return sorted_A;
    }
}

