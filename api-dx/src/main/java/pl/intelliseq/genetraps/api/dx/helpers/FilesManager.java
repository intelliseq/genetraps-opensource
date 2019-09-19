package pl.intelliseq.genetraps.api.dx.helpers;

import com.amazonaws.services.s3.AmazonS3;
import com.amazonaws.services.s3.AmazonS3ClientBuilder;
import com.amazonaws.services.s3.model.S3ObjectSummary;
import lombok.extern.log4j.Log4j2;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.core.env.Environment;
import pl.intelliseq.genetraps.api.dx.helpers.aws_manager.AWSApiProcessManager;

import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Created by intelliseq on 07/12/2017.
 */
@Log4j2
public class FilesManager {

    @Autowired
    private AWSApiProcessManager awsApiProcessManager;

    @Autowired
    private Environment env;

    private AtomicInteger counter = new AtomicInteger(1);

    final private AmazonS3 s3Client = AmazonS3ClientBuilder.defaultClient();

    public synchronized Integer mkdir() {
        var list = getNumericDirectories();
        if (list.contains(counter.get()) || list.size() != counter.get() - 1) {
            counter.set(getLowestFreeIndex(list));
        }
        awsApiProcessManager.runCreateSample(counter.get());
        return counter.getAndIncrement();
    }

    // AWS
    public List<Integer> getNumericDirectories() {
        String bucketName = env.getProperty("bucket.default");
        var listObjectsV2Result = s3Client.listObjectsV2(bucketName, String.format("%s/", env.getProperty("samples.folder")));
        List<S3ObjectSummary> sampleSummaries = listObjectsV2Result.getObjectSummaries();
        if(sampleSummaries.size() == 1) {
            return new LinkedList<>();
        }
        else {
            List<Integer> sampleIds = new ArrayList<>();
            for (S3ObjectSummary sampleSummary : sampleSummaries) {
                if(sampleSummary.getKey().matches("samples/([0-9]+)/")) {
                    Integer id = Integer.parseInt(sampleSummary.getKey().replaceFirst("samples/([0-9]+)/", "$1"));
                    sampleIds.add(id);
                }
            }
            Collections.sort(sampleIds);
            return sampleIds;
        }
    }

    public Integer getLowestFreeIndex(List<Integer> intList) {
        log.info("Getting lowest free index");
        //Jeśli nic nie ma to zwracamy 1
        if (intList == null || intList.size() == 0) {
            log.info("Lowest index: " + 1);
            return 1;
        } else if (intList.size() == intList.get(intList.size() - 1)) {
            //Jeśli ostatni element listy równa się jej wielkości, to oznacza, że nie ma dziur.
            log.info("Lowest index: " + (intList.size() + 1));
            return intList.size() + 1;
        }
        /**
         * W przeciwnym wypadku szukamy dziury
         */
        Integer index = IntStream
                .rangeClosed(1, intList.size())
                .filter(i -> !intList.contains(i))
                .sorted()
                .boxed()
                .collect(Collectors.toList())
                .get(0);
        log.info("Lowest index: " + index);
        return index;
    }
}
