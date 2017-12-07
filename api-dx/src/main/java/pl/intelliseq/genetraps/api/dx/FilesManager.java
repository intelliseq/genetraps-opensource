package pl.intelliseq.genetraps.api.dx;

import pl.intelliseq.genetraps.api.dx.parser.DxRunner;

import java.util.Objects;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Created by intelliseq on 07/12/2017.
 */
public class FilesManager {
    private AtomicInteger counter = new AtomicInteger(0);

    public synchronized Integer mkdir(){
        do{
            counter.incrementAndGet();
        }while(!Objects.equals(DxRunner.runCommand("dx mkdir " + counter.get()), ""));
        return counter.get();
    }

    public void test(){
        System.out.println("Hello World");
    }
}
