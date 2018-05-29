package pl.intelliseq.genetraps.api.dx.helpers;

import org.apache.log4j.Logger;
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

    private Logger log = Logger.getLogger(FilesManager.class);

    private AtomicInteger counter = new AtomicInteger(1);

    public synchronized Integer mkdir(){
        log.info("Counter: "+counter);
        String response = processManager.runCommand("dx ls --folders samples | grep -Eo "+counter);
        log.info("Response: "+response);
        if(!response.equals("")) {
            counter.set(getLowestFreeIndex());
            processManager.runMkdir(String.valueOf(counter.get()));
        }
        processManager.runMkdir(String.valueOf(counter.get()));
        return counter.getAndIncrement();
    }

    public Integer getAndIncrement(){
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
        log.info("Getting lowest free index");
        List<Integer> intList = getNumericDirectories();
        if(intList == null || intList.size() == 0){
            log.info("Lowest index: "+1);
            return 1;
        } else if(intList.size() == intList.get(intList.size()-1)){
            log.info("Lowest index: "+(intList.size()+1));
            return intList.size() +1;
        }
        /**
         * Bierzemy najniższy
         */
        Integer index = IntStream.rangeClosed(1, intList.size()).filter(i -> !intList.contains(i)).sorted().boxed().collect(Collectors.toList()).get(0);
        log.info("Lowest index: "+index);
        return index;
    }
}
