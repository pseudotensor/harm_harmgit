(custom-set-variables
  ;; custom-set-variables was added by Custom.
  ;; If you edit it by hand, you could mess it up, so be careful.
  ;; Your init file should contain only one such instance.
  ;; If there is more than one, they won't work right.
 '(inhibit-startup-screen t))
(custom-set-faces
  ;; custom-set-faces was added by Custom.
  ;; If you edit it by hand, you could mess it up, so be careful.
  ;; Your init file should contain only one such instance.
  ;; If there is more than one, they won't work right.
 )

(column-number-mode t)



(setq
  ;; don't show annoing startup msg
  inhibit-startup-message t
  ;; NO annoing backups
  make-backup-files nil
  ;; follow symlinks and don't ask
  vc-follow-symlinks t
  ;; see what you type
  echo-keystrokes 0.01
  ;; text scrolling
  scroll-conservatively 50
  scroll-preserve-screen-position 't
  ;; scroll-margin 10 ;; shifts too much and it's annoying to reject user choice of where want to put mouse or key
  ;; Insert space/slash after completion
  comint-completion-addsuffix t
  ;; number of chars in line
  ;fill-column 72
  ;; delete line in one stage
  kill-whole-line t
  ;; default mode
  default-major-mode 'text-mode
  ;; meaning are the same as the name :)
  delete-key-deletes-forward t
  ;; next-line don't add new lines
  next-line-add-newlines nil
  ;; Scroll by one line at a time
  scroll-step 1
  ;; don't add new lines when scrolling down
  next-line-add-newlines nil
  ;; make sure file ends with NEWLINE
  require-final-newline t
  ;; delete excess backup versions
  ;delete-old-versions t
  ;; setting the default tabulation
  default-tab-width 4
  ;; paste at cursor NOT at mouse pointer position
  mouse-yank-at-point t
  ;; apropos works better but slower
  apropos-do-all t
  ;; display time in the modeline
  display-time-24hr-format t
  display-time-day-and-date t
  ;; calendar customizing
  ;european-calendar-style t
  ;calendar-week-start-day 1
  ;; autosave every 512 keyboard inputs
  ;auto-save-interval 512
  ;; limit the number of newest versions
  ;kept-new-versions 5
  ;; limit the number of oldest versions
  ;kept-old-versions 5
  auto-save-list-file-prefix "~/.emacs.d/backups/save-"
  ;; don't beep
;  visible-bell t
  cursor-in-non-selected-windows nil
  ;; dired settings
  ;dired-recursive-copies t
  ;dired-recursive-deletes t
)

;; I hate tabs!
(setq-default indent-tabs-mode nil)

;; I hate the new default of cursor moving by visual lines.
(setq line-move-visual nil)

;;http://ergoemacs.org/emacs/emacs_editing_lisp.html

;; turn on highlight matching brackets when cursor is on one
(show-paren-mode 1)

(setq show-paren-style 'parenthesis) ; highlight just brackets
(setq show-paren-style 'expression) ; highlight entire bracket expression


;  From old .emacs file:

(setq auto-mode-alist
      (cons (cons "\\.m$" 'sm-mode)
          (cons (cons "\\.sav$" 'sm-mode) auto-mode-alist)))
(autoload 'sm-mode "~/sm/el/sm.el" nil t)



;; Visual feedback on selections
(setq-default transient-mark-mode t)
;; Enable wheelmouse support by default
(cond (window-system
       (mwheel-install)
))



;;; Matlab-mode setup:

;; Add local lisp folder to load-path
(setq load-path (append load-path (list "~/elisp")))

;; Set up matlab-mode to load on .m files
(autoload 'matlab-mode "matlab" "Enter MATLAB mode." t)
(setq auto-mode-alist (cons '("\\.m\\'" . matlab-mode) auto-mode-alist))
(autoload 'matlab-shell "matlab" "Interactive MATLAB mode." t)

;; Customization:
(setq matlab-indent-function t)	; if you want function bodies indented
(setq matlab-verify-on-save-flag nil) ; turn off auto-verify on save
(defun my-matlab-mode-hook ()
  (setq fill-column 76))		; where auto-fill should wrap
(add-hook 'matlab-mode-hook 'my-matlab-mode-hook)
(defun my-matlab-shell-mode-hook ()
  '())
(add-hook 'matlab-shell-mode-hook 'my-matlab-shell-mode-hook)

;; Turn off Matlab desktop
(setq matlab-shell-command-switches '("-nojvm"))

