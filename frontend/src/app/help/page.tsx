import React from 'react';

export default function HelpPage() {
  return (
    <div className="text-sm text-black leading-normal">
      <h2 className="text-xl font-bold border-b border-[#ccc] pb-1 mb-4 text-[#3A332D]">Help</h2>

      <p className="mb-6">
        This page provides detailed documentation on how to use the VaxElan 2.0 web server for multi-epitope vaccine design.
      </p>

      {/* Section 1 */}
      <h3 className="text-lg font-bold text-[#3A332D] mb-2">1. Input Data</h3>
      <hr className="mb-3 border-[#ccc]" />
      <p className="mb-3">
        VaxElan accepts pathogenic protein sequences in <strong>FASTA format</strong>. You may provide input in three ways:
      </p>
      <ul className="list-disc pl-6 mb-3 space-y-1">
        <li><strong>Upload FASTA File:</strong> Upload a .fasta, .fa, or .txt file containing one or more protein sequences.</li>
        <li><strong>Paste Sequence:</strong> Directly paste a FASTA-formatted sequence into the text area.</li>
        <li><strong>UniProt ID:</strong> Enter a valid UniProt accession number (e.g., P03135) and the server will automatically retrieve the sequence.</li>
      </ul>

      <p className="mb-2"><strong>FASTA Format Example:</strong></p>
      <pre className="bg-[#3A332D] text-[#F4F4EE] p-4 mb-4 text-xs font-mono border border-[#999] overflow-x-auto whitespace-pre">
{`>P03135 | Capsid protein VP1
MAPTVRSRSRSRSRSRSRSRSRSRSRSRSRSRSRSRSRS
SRSRSRSRSRSRSRSRSRSRSRSRSRSRSRSRSRSRSRS
...`}
      </pre>
      <p className="mb-6 text-stone-600">
        Ensure FASTA headers follow standard nomenclature. Multiple sequences are supported for batch assessment.
        Standard formats accepted: FASTA, UniProt ID.
      </p>

      {/* Section 2 */}
      <h3 className="text-lg font-bold text-[#3A332D] mb-2">2. Strategy Selection</h3>
      <hr className="mb-3 border-[#ccc]" />
      <p className="mb-3">
        After providing input, select one of the 6 available computational strategies. Each strategy targets a
        different research objective:
      </p>
      <table className="w-full border-collapse border border-[#999] mb-6 text-sm">
        <thead>
          <tr className="bg-[#3A332D] text-white">
            <th className="border border-[#999] px-3 py-2 text-left">Strategy</th>
            <th className="border border-[#999] px-3 py-2 text-left">Best For</th>
          </tr>
        </thead>
        <tbody>
          <tr><td className="border border-[#999] px-3 py-2 font-bold">Strategy 1: Full Pipeline</td><td className="border border-[#999] px-3 py-2">De novo discovery. Orchestrates MHC mapping, safety filtering, and structural assembly in a unified run.</td></tr>
          <tr className="bg-[#F4F4EE]"><td className="border border-[#999] px-3 py-2 font-bold">Strategy 3: Safety Priority</td><td className="border border-[#999] px-3 py-2">Already validated binders. Executes deep-tier toxicity and allergenicity assessments.</td></tr>
          <tr><td className="border border-[#999] px-3 py-2 font-bold">Strategy 6: Structural</td><td className="border border-[#999] px-3 py-2">3D modeling via ESMFold for structural stability and epitope exposure assessment.</td></tr>
        </tbody>
      </table>
      <p className="mb-6 text-stone-600">
        See the <a href="/algorithm" className="text-[#3A332D] underline font-bold">Algorithm</a> page for a complete description of all 6 strategies.
      </p>

      {/* Section 3 */}
      <h3 className="text-lg font-bold text-[#3A332D] mb-2">3. Interpreting Results</h3>
      <hr className="mb-3 border-[#ccc]" />
      <p className="mb-3">
        Once a job completes, results are displayed in the Results Viewer. Key output types include:
      </p>
      <ul className="list-disc pl-6 mb-4 space-y-1">
        <li><strong>Epitope Affinity Matrices:</strong> Quantifies sequence binding via IC50 (nM) or Percentile Rank. Focus on values below the 0.5% threshold for high-potency candidates.</li>
        <li><strong>Safety Diagnostic Flags:</strong> Heuristic alerts for toxicity and allergenicity. Any positive flag in ToxinPred assessments warrants immediate manual review or exclusion.</li>
        <li><strong>Structural Validation (PDB):</strong> In Strategy 6, 3D coordinate files (PDB) are generated via ESMFold. Visualize using PyMOL or ChimeraX to verify epitope exposure and linker stability.</li>
      </ul>
      <p className="mb-4">
        All result files (CSV, PDB) can be downloaded individually from the Results Viewer interface.
      </p>

      {/* Section 4 */}
      <h3 className="text-lg font-bold text-[#3A332D] mb-2">4. Execution Latency</h3>
      <hr className="mb-3 border-[#ccc]" />
      <p className="mb-4 bg-[#F4F4EE] p-3 border border-[#ccc]">
        <strong>Note:</strong> Deep structural modeling (Strategy 6) and aggregated filtered runs can exceed 20 minutes in latency. 
        The persistent job tracking system ensures data integrity throughout the asynchronous lifecycle. You may 
        close the browser and return to the Results Viewer using your Job ID.
      </p>

      {/* Reference */}
      <h3 className="text-lg font-bold text-[#3A332D] mb-2 border-t border-[#ccc] pt-4">Reference</h3>
      <p className="mb-6 bg-[#F4F4EE] p-3 border border-[#ccc]">
        If you are using VaxElan 2.0, please cite:<br/><br/>
        <em>Yukti, et al. &quot;VaxElan 2.0: A next-generation tool for prediction and designing of multi-epitope vaccines.&quot; (In Preparation, 2026)</em>
      </p>
    </div>
  );
}
