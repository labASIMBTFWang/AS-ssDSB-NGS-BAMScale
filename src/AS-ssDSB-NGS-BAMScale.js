// @ts-check

// conda install -c bioconda bamscale
// conda install -c bioconda ucsc-bigwigtobedgraph
// conda install -c bioconda ucsc-bigwigtowig

const path = require("path");
const fs = require("fs");
const child_process = require("child_process");

const MAX_THREAD = 16;

/**
 * @template T
 * @param {string[]} argv
 * @param {RegExp} regexp
 * @param {(arg: string, args: string[]) => T} trx
 * @param {T} [default_val]
 * @returns {T}
 * @example
 * get(/--aaa=(\d+)/, (arg, args) => Number(args[0]))
 */
function getArg(argv, regexp, trx, default_val) {
	const aa = argv.map(s => s.match(regexp)).filter(a => a)[0];
	if (aa) {
		const i = argv.indexOf(aa[0]);
		if (i >= 0) {
			const args = argv[i];

			const r = trx?.(aa[0], aa.slice(1));

			argv.splice(i, 1);

			return r;
		}
	}
	return /** @type {any} */(default_val);
}

/**
 * @param {import("stream").Readable} stream
 * @returns {Promise<Buffer>}
 */
async function streamToBuffer(stream) {
	if (!stream) {
		throw new TypeError();
	}
	return new Promise((resolve, reject) => {
		const chunks = [];
		stream.on("error", (err) => reject(err));
		stream.on("data", (chunk) => chunks.push(Buffer.from(chunk)));
		stream.on("end", () => resolve(Buffer.concat(chunks)));
	});
}

/**
 * @param {import("stream").Readable} stream
 * @param {any} [encoding]
 * @returns {Promise<string>}
 */
async function streamToString(stream, encoding) {
	return (await streamToBuffer(stream)).toString(encoding);
}

class ProcessUtiil {
	/**
	 * exitPromise
	 * @type {Promise<child_process.ChildProcess>}
	 */
	promise;
	
	/**
	 * do not await proc.promise // process no write stdio if exit
	 */
	async stdoutToString() { return ""; }
	/**
	 * do not await proc.promise // process no write stdio if exit
	 */
	async stderrToString() { return ""; }
	
	/**
	 * @param {string} filename
	 * @returns {fs.WriteStream|null}
	 */
	stdoutToFile(filename) { return null; }
	/**
	 * @param {string} filename
	 * @returns {fs.WriteStream|null}
	 */
	stderrToFile(filename) { return null; }
	/**
	 * @param {string} out_filename
	 * @param {string} err_filename
	 * @returns {(fs.WriteStream|null)[]}
	 */
	stdioToFile(out_filename, err_filename) { return [null, null, null]; }

	/**
	 * @param {import("child_process").ChildProcess} proc
	 * @param {Promise<child_process.ChildProcess>} exit_promise
	 */
	static assign(proc, exit_promise) {
		const stdoutToString = async () => proc.stdout ? await streamToString(proc.stdout) : ""; // do not await proc.promise // process write stdio when live
		const stderrToString = async () => proc.stderr ? await streamToString(proc.stderr) : ""; // do not await proc.promise // process write stdio when live
		/** @param {string} filename */
		const stdoutToFile = (filename) => proc.stdout ? proc.stdout.pipe(fs.createWriteStream(filename)) : null;
		/** @param {string} filename */
		const stderrToFile = (filename) => proc.stderr ? proc.stderr.pipe(fs.createWriteStream(filename)) : null;
		/**
		 * @param {string} out_filename
		 * @param {string} err_filename
		 */
		const stdioToFile = (out_filename, err_filename) => {
			return [
				null,
				stdoutToFile(out_filename),
				stderrToFile(err_filename),
			];
		};

		/**
		 * @type {ProcessUtiil}
		 */
		const uu = {
			/** exit_promise */
			promise: exit_promise,
			stdoutToString: stdoutToString,
			stderrToString: stderrToString,
			stdoutToFile: stdoutToFile,
			stderrToFile: stderrToFile,
			stdioToFile: stdioToFile,
			// then: promise.then,
			// catch: promise.catch,
		};
		
		/** @type {child_process.ChildProcess & ProcessUtiil} */
		const v = Object.assign(proc, uu)
		return v;
	}
}

execAsync.log = "log";
/**
 * @param {string} cmd
 * @param {child_process.ExecOptions} [options]
 * returns {child_process.ChildProcess & { promise: Promise<child_process.ChildProcess>; then: typeof Promise.prototype.then; catch: typeof Promise.prototype.catch; stdoutToString: () => Promise<string>; stderrToString: () => Promise<string>; }}
 * @returns {child_process.ChildProcess & ProcessUtiil}
 */
function execAsync(cmd, options) {
	console[execAsync.log]("exec:", cmd);

	const proc = options ? child_process.exec(cmd, options) : child_process.exec(cmd);

	const forward_err = new Error(cmd);// forward error
	
	const promise = new Promise((resolve, reject) => {
		proc.on("error", err => {
			// @ts-ignore
			forward_err.src = err.stack;
			reject(forward_err);
		});
		proc.on("exit", (code, signal) => resolve(proc));
	});

	return ProcessUtiil.assign(proc, promise);
}

/**
 * @param {string} dir
 * @param {string} prefix
 * @param {string} ref_fa
 * @param {string} fq1
 * @param {string} fq2
 * @param {number} max_thread
 */
async function run_bwa_mem(dir, prefix, ref_fa, fq1, fq2, max_thread, { run = false }) {
	const logDir = path.resolve("./", dir, "logs");
	const outFile = path.resolve(dir, `${prefix}.sort.bam`);

	const cmd = `bwa mem ${ref_fa} ${fq1} ${fq2} | samtools sort -@ ${max_thread} -l 9 -o ${outFile}`;
	
	if (run) {
		if (!fs.existsSync(logDir)) {
			fs.mkdirSync(logDir);
		}

		const proc = execAsync(cmd);
		proc.stdioToFile(`${logDir}/${prefix}.stdout.txt`, `${logDir}/${prefix}.stderr.txt`);
		await proc.promise;

		await execAsync(`samtools index ${outFile}`).promise;
	}
	else {
		console.error(cmd);
	}

	return outFile;
}

/**
 * @param {string} sampleName
 * @param {string} scale_bw
 * @param {{ run: boolean }} _
 * @returns {Promise<string>} out_file
 * bedGraph format: ["chr", "start", "end", "real"]
 */
async function run_bigWigToBed(sampleName, scale_bw, { run = false }) {
	// const out_file = `${sampleName}.sort.bam.scaled.bed`;
	const out_file = scale_bw.replace(/\.bw$/, ".bed");

	const bedGraph = scale_bw.replace(/\.bw$/, ".bedGraph"); // IGV
	const wig = scale_bw.replace(/\.bw$/, ".wig");

	const cmd_bigWigToBedGraph = `bigWigToBedGraph ${scale_bw} ${bedGraph}`
	const cmd_bigWigToWig = `bigWigToWig ${scale_bw} ${wig}`
	const cmd_wig2bed = `wig2bed < ${wig} > ${out_file}`

	if (run) {
		await execAsync(cmd_bigWigToBedGraph).promise;
		
		await execAsync(cmd_bigWigToWig).promise;
		await execAsync(cmd_wig2bed).promise;
	}
	else {
		console.error(cmd_bigWigToBedGraph);
		console.error(cmd_bigWigToWig);
		console.error(cmd_wig2bed);
	}
	
	return out_file;
}

/**
 * 
 * @param {string} sample_name
 * @param {string} bam_path
 * @param {string} operation default: scaled
 * @param {number} binSize default: 10
 * @param {{ run: boolean }} _
 * @returns {Promise<string>} output file path ${input_file_dir}/${BAMscale args...}/${input_file_name}.scaled.bw
 */
async function run_BAMscale(sample_name, bam_path, operation, binSize, { run = false }) {
	const input_file_dir = path.dirname(bam_path);
	const input_file_name = path.basename(bam_path);
	const outdir = `${input_file_dir}/${sample_name}_bin${binSize}_${operation}`;

	// const bai = `${bam_path}.bai`;
	const cmd_bai = `samtools index ${bam_path}`;
	const cmd_scale = `BAMscale scale -t ${MAX_THREAD} --operation ${operation} --binsize ${binSize} --outdir ${outdir} --bam ${bam_path}`;

	if (run) {
		await execAsync(cmd_bai).promise;
		await execAsync(cmd_scale).promise;
	}
	else {
		console.error(cmd_bai);
		console.error(cmd_scale);
	}

	return `${outdir}/${input_file_name}.scaled.bw`;
}

async function main() {
	const path_to_dataset = getArg(process.argv, /--samples=(.+)/, (arg, args) => args[0]);
	const binSize = getArg(process.argv, /--bin=(.+)/, (arg, args) => Number(args[0]), 10);
	const run = !getArg(process.argv, /--dry-run/, (arg, args) => true);
	
	if (path_to_dataset && binSize) {
		const dataset = JSON.parse(fs.readFileSync(path_to_dataset).toString());

		const ppp = dataset.map(async ref_param => {
			const pp = ref_param.samples.map(async sample_param => {
				const sample_name = [ref_param.reference_name, sample_param.sample_name].join("_");
				const [fq1, fq2] = sample_param.reads;
				
				const bam = await run_bwa_mem("./", sample_name, ref_param.reference_genome_fasta, fq1, fq2, MAX_THREAD, { run: run, });
				
				const scale_bw = await run_BAMscale(sample_name, bam, "scaled", binSize, { run: run, });
		
				return await run_bigWigToBed(sample_name, scale_bw, { run: run, });
			});
			return await Promise.all(pp);
		});
		console.log((await Promise.all(ppp)).flat(1).join("\n"));
	}
	else {
		console.log({
			path_to_dataset,
			binSize,
			run,
		});
		console.log(`example: node ssDNA_peak_pipeline_public.js --bin=10 --samples=example_samples.json`);
		console.log(`print cmd: node ssDNA_peak_pipeline_public.js --dry-run --bin=10 --samples=example_samples.json`);
	}
}

main();

