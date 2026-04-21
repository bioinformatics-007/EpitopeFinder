import React from 'react';

export default function AlgorithmPage() {
  return (
    <div className="text-sm text-black leading-normal">
      <h2 className="text-xl font-bold border-b border-[#ccc] pb-1 mb-4 text-[#3A332D]">Algorithm</h2>

      <p className="mb-4">
        VaxElan 2.0 implements <strong>6 specialized computational strategies</strong> for multi-epitope vaccine design. 
        Each strategy orchestrates a distinct combination of prediction, filtering, and assembly tools to address 
        different research objectives. The table below summarizes all available strategies.
      </p>

      <table className="w-full border-collapse border border-[#999] mb-6 text-sm">
        <thead>
          <tr className="bg-[#3A332D] text-white">
            <th className="border border-[#999] px-3 py-2 text-left w-12">#</th>
            <th className="border border-[#999] px-3 py-2 text-left w-48">Strategy Name</th>
            <th className="border border-[#999] px-3 py-2 text-left">Description</th>
            <th className="border border-[#999] px-3 py-2 text-left w-56">Tools Used</th>
          </tr>
        </thead>
        <tbody>
          <tr>
            <td className="border border-[#999] px-3 py-2 font-bold">1</td>
            <td className="border border-[#999] px-3 py-2 font-bold">Multi-epitope Pipeline</td>
            <td className="border border-[#999] px-3 py-2">The standard end-to-end pipeline for designing multivalent vaccines. Includes B-cell epitope prediction, CTL (MHC-I) and HTL (MHC-II) prediction, followed by population coverage analysis and assembly.</td>
            <td className="border border-[#999] px-3 py-2">IEDB MHC-I/II, NetCTL, B-cell (Kolaskar &amp; Tongaonkar)</td>
          </tr>
          <tr className="bg-[#F4F4EE]">
            <td className="border border-[#999] px-3 py-2 font-bold">2</td>
            <td className="border border-[#999] px-3 py-2 font-bold">Targeted Assessment</td>
            <td className="border border-[#999] px-3 py-2">A modular strategy allowing for the execution of specific tools on a per-pathogen basis. Ideal for researchers who have already identified candidate proteins and need precise epitope mapping.</td>
            <td className="border border-[#999] px-3 py-2">NetCTL-Pan, MHC-I Binding, MHC-II Binding</td>
          </tr>
          <tr>
            <td className="border border-[#999] px-3 py-2 font-bold">3</td>
            <td className="border border-[#999] px-3 py-2 font-bold">Safety &amp; Host Compatibility</td>
            <td className="border border-[#999] px-3 py-2">Focuses on the safety profile of vaccine candidates. Filters out toxic, allergenic, or human-homologous peptides to minimize the risk of autoimmune responses.</td>
            <td className="border border-[#999] px-3 py-2">ToxinPred2, DeepTMHMM, SignalP, Human Homology Check</td>
          </tr>
          <tr className="bg-[#F4F4EE]">
            <td className="border border-[#999] px-3 py-2 font-bold">4</td>
            <td className="border border-[#999] px-3 py-2 font-bold">Integrated Prediction</td>
            <td className="border border-[#999] px-3 py-2">Combines multiple predictive models to identify epitopes that are simultaneously potent and safe. Runs Strategies 1, 2, and 3 in parallel and aggregates the best candidates.</td>
            <td className="border border-[#999] px-3 py-2">All Top-tier tools, Aggregated Ranking</td>
          </tr>
          <tr>
            <td className="border border-[#999] px-3 py-2 font-bold">5</td>
            <td className="border border-[#999] px-3 py-2 font-bold">Sequential Epitope Filtering</td>
            <td className="border border-[#999] px-3 py-2">Uses a sequential approach: identify epitopes, then filter for allergenicity, then toxicity, then antigenicity. Ensures every resulting sequence meets all selection criteria.</td>
            <td className="border border-[#999] px-3 py-2">AlgPred, IAPred, ToxinPred3</td>
          </tr>
          <tr className="bg-[#F4F4EE]">
            <td className="border border-[#999] px-3 py-2 font-bold">6</td>
            <td className="border border-[#999] px-3 py-2 font-bold">3D Structural Validation</td>
            <td className="border border-[#999] px-3 py-2">The final step in vaccine design. Uses ESMFold to predict the 3D structure of the assembled vaccine construct, allowing for assessment of structural stability and epitope exposure.</td>
            <td className="border border-[#999] px-3 py-2">ESMFold, SASA Filter</td>
          </tr>
        </tbody>
      </table>

      <h3 className="text-lg font-bold text-[#3A332D] mb-2 border-t border-[#ccc] pt-4">Methodology</h3>
      <p className="mb-4">
        Our platform integrates the latest benchmarks from the IEDB and leading computational models to ensure 
        that every predicted epitope is validated across binding affinity, sequence conservation, and host compatibility. 
        Heuristic modeling, structural validation, and clinical safety checks are applied at each stage.
      </p>

      <h3 className="text-lg font-bold text-[#3A332D] mb-2 border-t border-[#ccc] pt-4">Reference</h3>
      <p className="mb-6 bg-[#F4F4EE] p-3 border border-[#ccc]">
        If you are using VaxElan 2.0, please cite:<br/><br/>
        <em>Yukti, et al. &quot;VaxElan 2.0: A next-generation tool for prediction and designing of multi-epitope vaccines.&quot; (In Preparation, 2026)</em>
      </p>
    </div>
  );
}
