package pl.intelliseq.genetraps.api.dx;

import pl.intelliseq.genetraps.api.dx.parser.DxRunner;

import java.util.Arrays;
import java.util.List;
import java.util.Objects;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * Created by intelliseq on 07/12/2017.
 */
public class FilesManager {
    private AtomicInteger counter = new AtomicInteger(1);

    public synchronized Integer mkdir(){
        String response = DxRunner.runCommand("dx mkdir " + counter.get());
        if(!response.equals("")) {
            counter.set(getLowestFreeIndex());
            DxRunner.runCommand("dx mkdir " + counter.get());
        }
        return counter.getAndIncrement();
    }

    public void resetCounter(){
        counter.set(1);
    }

    public List<Integer> getNumericDirectories(){
        String command = "dx ls --folders | grep -Eo \"^([0-9]+)\"";
        String response = DxRunner.runCommand(command);
        if(response.equals("")){
            return null;
        }
        List<Integer> list = Arrays.stream(response.split("\n")).map(Integer::valueOf).sorted().collect(Collectors.toList());
        return list;
    }

    public Integer getLowestFreeIndex(){
        List<Integer> intList = getNumericDirectories();
        if(intList == null){
            return 0;
        } else if(intList.size() == intList.get(intList.size()-1)){
            return intList.size() +1;
        }
        return IntStream.rangeClosed(1, intList.size()).filter(i -> !intList.contains(i)).sorted().boxed().collect(Collectors.toList()).get(0);
    }

    public void test(){
        System.out.println("Hello World");
    }
}
