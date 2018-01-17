package pl.intelliseq.genetraps.api.dx.enums;

public enum JobStates {
    IDLE("idle"),
    RUNNABLE("runnable"),
    RUNNING("running"),
    DONE("done"),
    WAITING_ON_INPUT("waiting_for_input"),
    WAITING_ON_OUTPUT("waiting_on_output");

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
        throw new IllegalArgumentException();
    }
}
