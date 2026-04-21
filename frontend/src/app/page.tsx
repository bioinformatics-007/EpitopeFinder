import React from 'react';
import Link from 'next/link';

export default function LandingPage() {
  return (
    <div className="text-sm text-black leading-normal">
      <h2 className="text-xl font-bold border-b border-[#ccc] pb-1 mb-4 text-[#3A332D]">Welcome to VaxElan 2.0 Web Server</h2>
      
      <p className="mb-4">
        <strong>VaxElan</strong> is an <em>in silico</em> method developed to predict and design multi-epitope vaccine constructs from pathogenic protein sequences. This server integrates robust machine learning classifiers, IEDB stabilization metrics, and structural assessment protocols to prioritize optimal vaccine candidates.
      </p>

      <p className="mb-4">
        Our comprehensive pipeline natively supports both B-cell and T-cell (MHC-I, MHC-II) epitope prediction, filtering candidates based on toxicity, allergenicity, and antigenicity. Finally, the server links prioritized epitopes using validated structural linkers and predicts the 3D conformation of the final multivalent vaccine using ESMFold integration.
      </p>

      <h3 className="text-lg font-bold text-[#3A332D] mt-6 mb-2">Major Features</h3>
      <ul className="list-disc pl-6 mb-6 space-y-1">
        <li><strong>Design Vaccine:</strong> Submit pathogenic FASTA sequences for end-processing.</li>
        <li><strong>Toxicity & Allergenicity:</strong> Utilize strict non-toxic and non-allergenic filtering (inspired by ToxinPred).</li>
        <li><strong>Structural Modeling:</strong> High-fidelity 3D modeling using ESMFold algorithms.</li>
        <li><strong>Comprehensive Results:</strong> Download all node CSVs or view analytical insights natively.</li>
      </ul>

      <h3 className="text-lg font-bold text-[#3A332D] mb-2 border-t border-[#ccc] pt-4">Reference</h3>
      <p className="mb-6 bg-[#F4F4EE] p-3 border border-[#ccc]">
        If you are using VaxElan 2.0, please cite:<br/><br/>
        <em>Yukti, et al. "VaxElan 2.0: A next-generation tool for prediction and designing of multi-epitope vaccines." (In Preparation, 2026)</em>
      </p>

      <div className="flex gap-4">
        <Link href="/submit" className="border border-[#999] bg-[#E0E0E0] px-4 py-1 text-black hover:bg-[#ccc] no-underline">
          Click Here to Start Design
        </Link>
        <Link href="/help" className="border border-[#999] bg-[#E0E0E0] px-4 py-1 text-black hover:bg-[#ccc] no-underline">
          View Protocol Help
        </Link>
      </div>
    </div>
  );
}
