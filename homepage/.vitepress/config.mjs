import { defineConfig } from "vitepress";

// https://vitepress.dev/reference/site-config
export default defineConfig({
  lang: "en-US",
  title: "Lassa GPC DMS",
  description:
    "Interactive figures and detailed results for deep mutational scanning of the GPC from the lineage IV Lassa virus Josiah strain.",
  base: "/LASV_Josiah_GP_DMS",
  appearance: false,
  themeConfig: {
    // https://vitepress.dev/reference/default-theme-config
    nav: [
      { text: "Home", link: "/" },
      { text: 'Paper', link: "https://www.cell.com/immunity/fulltext/S1074-7613(24)00319-4" },
      { text: 'Cell Entry', link: '/cell_entry' },
      { text: 'Antibody Escape', link: '/antibody_escape' },
      { text: "Appendix", link: "/appendix", target: "_self" },
    ],
    socialLinks: [{ icon: "github", link: "https://github.com/dms-vep/LASV_Josiah_GP_DMS.git" }],
    footer: {
      message: 'Study by <a href="https://www.cell.com/immunity/fulltext/S1074-7613(24)00319-4">Carr, Crawford, et al (2024)</a>',
    },
  },
});
