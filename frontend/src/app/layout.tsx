import type { Metadata } from 'next';
import './globals.css';
import { HeaderBanner } from '@/components/layout/HeaderBanner';
import { HorizontalMenu } from '@/components/layout/HorizontalMenu';
import { Footer } from '@/components/layout/Footer';

export const metadata: Metadata = {
  title: 'VaxElan — Computational Vaccine Design Platform',
  description:
    'Design multi-epitope vaccines with automated B-cell, MHC-I, and MHC-II epitope prediction, protein prioritization, and 3D structural validation.',
  keywords: [
    'vaccine design',
    'epitope prediction',
    'bioinformatics',
    'multi-epitope vaccine',
    'MHC prediction',
    'computational immunology',
  ],
};

export default function RootLayout({
  children,
}: {
  children: React.ReactNode;
}) {
  return (
    <html lang="en">
      <body className="font-sans bg-[#F4F4EE] text-black m-0">
        <HeaderBanner />
        <HorizontalMenu />
        <div className="w-full max-w-[1000px] mx-auto bg-white border-x border-b border-[#999] min-h-[500px]">
          <main className="p-4">
             {children}
          </main>
        </div>
        <Footer />
      </body>
    </html>
  );
}
