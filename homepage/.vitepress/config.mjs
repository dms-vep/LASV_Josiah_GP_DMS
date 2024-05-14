import { defineConfig } from "vitepress";

// https://vitepress.dev/reference/site-config
export default defineConfig({
  lang: "en-US",
  title: "Deep mutational scanning of Lassa virus Josiah strain glycoprotein complex (GPC)",
  description:
    "Interactive figures and detailed results for deep mutational scanning of the GPC from the lineage IV Lassa virus Josiah strain.",
  base: "/LASV_Josiah_GP_DMS",
  appearance: false,
  themeConfig: {
    // https://vitepress.dev/reference/default-theme-config
    nav: [
      { text: "Home", link: "/" },
      { text: "Appendix", link: "/appendix", target: "_self" },
    ],
    socialLinks: [{ icon: "github", link: "https://github.com/dms-vep/LASV_Josiah_GP_DMS.git" }],
    footer: {
      message: "Study by Caleb Carr, Kate Crawford, Jesse Bloom, et al",
    },
  },
});
