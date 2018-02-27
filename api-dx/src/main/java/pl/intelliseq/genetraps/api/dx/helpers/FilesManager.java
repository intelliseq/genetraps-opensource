package pl.intelliseq.genetraps.api.dx.helpers;

import org.springframework.beans.factory.annotation.Autowired;

import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Created by intelliseq on 07/12/2017.
 */
public class FilesManager {
    @Autowired
    private ProcessManager processManager;

    private AtomicInteger counter = new AtomicInteger(1);

    public synchronized Integer mkdir(){
        String response = processManager.runMkdir(String.valueOf(counter.get()));
        if(!response.equals("")) {
            counter.set(getLowestFreeIndex());
            processManager.runMkdir(String.valueOf(counter.get()));
        }
        return counter.getAndIncrement();
    }

    public void resetCounter(){
        counter.set(1);
    }

    public List<Integer> getNumericDirectories(){
        String command = "dx ls --folders samples | grep -Eo \"^([0-9]+)\"";
        String response = processManager.runCommand(command);
        if(response.equals("")){
            return new LinkedList<>();
        }
        return Arrays.stream(response.split("\n")).map(Integer::valueOf).sorted().collect(Collectors.toList());
    }

    public Integer getLowestFreeIndex(){
        List<Integer> intList = getNumericDirectories();
        if(intList == null || intList.size() == 0){
            return 0;
        } else if(intList.size() == intList.get(intList.size()-1)){
            return intList.size() +1;
        }
        /**
         * Bierzemy najniższy
         */
        return IntStream.rangeClosed(1, intList.size()).filter(i -> !intList.contains(i)).sorted().boxed().collect(Collectors.toList()).get(0);
    }
}
