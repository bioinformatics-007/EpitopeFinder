import React from 'react';

export default function ContactPage() {
  return (
    <div className="text-sm text-black leading-normal">
      <h2 className="text-xl font-bold border-b border-[#ccc] pb-1 mb-4 text-[#3A332D]">Contact</h2>

      <p className="mb-4">
        For technical queries, collaboration requests, or bug reports related to EpitopePred, please contact the development team using the information below.
      </p>

      <table className="w-full border-collapse border border-[#999] mb-6 text-sm">
        <tbody>
          <tr className="bg-[#3A332D] text-white">
            <th className="border border-[#999] px-3 py-2 text-left" colSpan={2}>Development Team</th>
          </tr>
          <tr>
            <td className="border border-[#999] px-3 py-2 font-bold w-40">Team</td>
            <td className="border border-[#999] px-3 py-2">EpitopePred Development Team — Bioinformatics &amp; Software Engineering</td>
          </tr>
          <tr className="bg-[#F4F4EE]">
            <td className="border border-[#999] px-3 py-2 font-bold">Description</td>
            <td className="border border-[#999] px-3 py-2">A multidisciplinary team dedicated to accelerating vaccine discovery through automated computational pipelines and advanced biological heuristics.</td>
          </tr>
          <tr>
            <td className="border border-[#999] px-3 py-2 font-bold">Email</td>
            <td className="border border-[#999] px-3 py-2"><a href="mailto:bioinformatics601@gmail.com" className="text-[#3A332D] underline">bioinformatics601@gmail.com</a></td>
          </tr>
          <tr className="bg-[#F4F4EE]">
            <td className="border border-[#999] px-3 py-2 font-bold">Source Code</td>
            <td className="border border-[#999] px-3 py-2"><a href="https://github.com/bioinformatics-007/EpitopePred" target="_blank" rel="noreferrer" className="text-[#3A332D] underline">github.com/bioinformatics-007/EpitopePred</a></td>
          </tr>
        </tbody>
      </table>

      <h3 className="text-lg font-bold text-[#3A332D] mb-2 border-t border-[#ccc] pt-4">Related Resources</h3>
      <p className="mb-4">
        EpitopePred is built upon the following global biological databases and prediction servers:
      </p>
      <ul className="list-disc pl-6 mb-6 space-y-1">
        <li><a href="https://www.iedb.org/" target="_blank" rel="noreferrer" className="text-[#3A332D] underline">IEDB — Immune Epitope Database</a></li>
        <li><a href="https://www.uniprot.org/" target="_blank" rel="noreferrer" className="text-[#3A332D] underline">UniProt — Universal Protein Knowledgebase</a></li>
        <li><a href="https://webs.iiitd.edu.in/raghava/toxinpred/" target="_blank" rel="noreferrer" className="text-[#3A332D] underline">ToxinPred — Prediction of Toxic Peptides</a></li>
        <li><a href="https://webs.iiitd.edu.in/raghava/algpred/" target="_blank" rel="noreferrer" className="text-[#3A332D] underline">AlgPred — Allergenicity Prediction</a></li>
      </ul>

      <h3 className="text-lg font-bold text-[#3A332D] mb-2 border-t border-[#ccc] pt-4">Reference</h3>
      <p className="mb-6 bg-[#F4F4EE] p-3 border border-[#ccc]">
        If you are using EpitopePred, please cite:<br/><br/>
        <em>&quot;EpitopePred: A next-generation tool for prediction and designing of multi-epitope vaccines.&quot;</em>
      </p>
    </div>
  );
}
