import React from 'react';
import Link from 'next/link';

export function HeaderBanner() {
  return (
    <header className="w-full bg-[#3A332D] shrink-0">
      <div className="mx-auto max-w-[1000px] py-6 px-4">
        <Link href="/" className="block no-underline">
          <h1 className="text-4xl font-bold text-white m-0 mb-1" style={{ fontFamily: 'Georgia, "Times New Roman", serif' }}>
            VaxElan 2.0
          </h1>
          <h2 className="text-lg font-normal text-stone-300 m-0" style={{ fontFamily: 'Georgia, "Times New Roman", serif' }}>
            Designing and Prediction of Multi-Epitope Vaccines
          </h2>
        </Link>
      </div>
    </header>
  );
}
