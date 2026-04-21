/** @type {import('tailwindcss').Config} */
module.exports = {
  content: [
    './src/app/**/*.{js,ts,jsx,tsx,mdx}',
    './src/components/**/*.{js,ts,jsx,tsx,mdx}',
  ],
  darkMode: 'class',
  theme: {
    extend: {
      colors: {
        amber: {
          50: '#F4F4EE',
          100: '#EBEBE2',
          200: '#D5D6C7',
          300: '#BDBFA7',
          400: '#A4AA04',
          500: '#8A8F03',
          600: '#635B53',
          700: '#4A433A',
          800: '#3A332D',
          900: '#26211D',
          950: '#181412',
        },
      },
      fontFamily: {
        sans: ['Arial', 'Helvetica', 'sans-serif'],
        mono: ['"Courier New"', 'Courier', 'monospace'],
      },
      maxWidth: {
        content: '1000px',
      },
    },
  },
  plugins: [],
};
