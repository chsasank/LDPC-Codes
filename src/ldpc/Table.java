
package ldpc;

import java.util.*;
/**
 * Basically a hash map with multiple keys
 * @author Sasank Chilamkurthy
 */
public class Table<V> {
    private Map<Integer, Map<Integer, V>> map;
    
    public Table(){
        map = new HashMap<>();
    }
    
    public void put(int i, int j, V v){
        if(map.containsKey(i)) map.get(i).put(j, v);
        else{
            Map x = new HashMap<>();
            x.put(j,v);
            map.put(i, x);
        }
    }
    
    public V get(int i, int j){
        return map.get(i).get(j);
    }
}
