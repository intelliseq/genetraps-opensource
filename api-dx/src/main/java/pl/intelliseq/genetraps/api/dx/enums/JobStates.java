package pl.intelliseq.genetraps.api.dx.enums;

import org.apache.log4j.Logger;

public enum JobStates {
    IDLE("idle"),
    RUNNABLE("runnable"),
    RUNNING("running"),
    DONE("done"),
    WAITING_ON_INPUT("waiting_for_input"),
    WAITING_ON_OUTPUT("waiting_on_output");

    private static Logger log = Logger.getLogger(JobStates.class);

    private String name;

    JobStates(String name){
        this.name = name;
    }

    @Override
    public String toString(){
        return name;
    }

    public static JobStates getEnum(String name){
        for(JobStates jobStates:values()){
            if (jobStates.name.equals(name)){
                return jobStates;
            }
        }
        log.error("IllegalArgumentException: "+name);
        throw new IllegalArgumentException();
    }
}
