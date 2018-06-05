package pl.intelliseq.genetraps.api.dx.helpers;

import com.dnanexus.DXContainer;
import com.dnanexus.exceptions.ResourceNotFoundException;
import org.apache.log4j.Logger;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.core.env.Environment;

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
    private DxApiProcessManager processManager;

    @Autowired
    private Environment env;

    private Logger log = Logger.getLogger(FilesManager.class);

    private AtomicInteger counter = new AtomicInteger(1);

    public synchronized Integer mkdir(){
        if(getNumericDirectories().contains(counter.get())) {
            counter.set(getLowestFreeIndex());
        }
        processManager.runMkDir(counter.get());
        return counter.getAndIncrement();
    }

    public Integer getAndIncrement(){
        return counter.getAndIncrement();
    }

    public void resetCounter(){
        counter.set(1);
    }

    public List<Integer> getNumericDirectories(){
        try {
            return DXContainer.getInstance(env.getProperty("dx-project"))
                    .listFolder("/samples")
                    .getSubfolders()
                    .stream()
                    .map(s -> Integer.parseInt(s.substring(9)))
                    .sorted()
                    .collect(Collectors.toList());
        }catch (ResourceNotFoundException e){
            return new LinkedList<>();
        }
    }

    public Integer getLowestFreeIndex(){
        log.info("Getting lowest free index");
        List<Integer> intList = getNumericDirectories();
        System.out.println(intList);
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
