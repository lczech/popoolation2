import java.util.ArrayList;
import java.util.concurrent.Callable;
import java.util.concurrent.*;
public class BatchProcessor 
{
	private final String[] batch;
	private final ExecutorService executor;
	private final MpileupParser parser;
	
	public BatchProcessor(String[] batch,MpileupParser parser, ExecutorService executor){
		this.batch=batch;
		this.executor=executor;
		this.parser=parser;
	}
	
	public String[] processBatch(){
		
		SyncBatchCollector batchCol=new SyncBatchCollector(this.batch.length);
		
		ArrayList<Callable<Object>> call=new ArrayList<Callable<Object>>();
		for(int i=0; i < this.batch.length; i++)
		{
			call.add(Executors.callable(new ProcessSingleEntry(this.batch[i],i,batchCol, parser)));
		}
		

		try
		{	
			// Run them all!
			 executor.invokeAll(call);	
		}
		catch(InterruptedException e)
		{
			e.printStackTrace();
			System.exit(0);
		}
		String[] batch=batchCol.getBatch();
		for(String s:batch)
		{
			if(s==null)
			{
					executor.shutdownNow();
					throw new IllegalArgumentException("An error occured during parsing of the mpileup;\n"+
							"Please make sure you are a.) using the correct encoding and b.) using the correct mpileup format");
			}
			}
			
		return batch;
	}

}


/**
 * Class for collecting the processed output
 * @author robertkofler
 *
 */
class SyncBatchCollector
{
	private String[] batch;
	public SyncBatchCollector(int blockSize)
	{
		batch=new String[blockSize];
	}
	
	public synchronized void addEntry(String entry, int index)
	{
		this.batch[index]=entry;
	}
	public synchronized String[] getBatch()
	{
		return this.batch;
	}
}

/**
 * Perform a single conversion event; Parse the mpileup and convert it into a sync
 * @author robertkofler
 *
 */
class ProcessSingleEntry implements Runnable
{
	private String toProcess;
	private int index;
	private SyncBatchCollector batchCol;
	private MpileupParser parser;
	public ProcessSingleEntry(String toProcess,int index, SyncBatchCollector batchCol,MpileupParser parser)
	{
		this.toProcess=toProcess;
		this.index=index;
		this.batchCol=batchCol;
		this.parser =parser;
	}
	
	@Override
	public void run()
	{
		SyncLine sl=parser.parseLine(this.toProcess);
		String formated=SyncFormater.formatSync(sl);
		batchCol.addEntry(formated, this.index);
	}
}
